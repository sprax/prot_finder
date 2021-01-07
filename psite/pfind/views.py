# @file: pfind/views.py
'''
Django views / form-page controllers.
NOTE: Function whose names start with a verb are procedural.
      Function whose names start with a noun or adjective and noun
        generally return that kind of object.

ncbi_name: NCBI Protein name as primary key.
found_idx: Smallest offset index found (so far) for a search term per Protein.
    - Subject to change because if 'acgtacgt' is found at offset X,
        it implies that the substring 'tacg' was found at offset X + 3;
        but 'tacg' might also be found at a smaller offset later, by a
        different search.
'''
from asgiref.sync import sync_to_async
from django.contrib.auth.models import User
# from django.core import serializers
from django.db import IntegrityError
from django.db.models import Count, Value
from django.db.models.functions import StrIndex
# from django.forms.models import model_to_dict
from django.http import Http404
from django.http import HttpResponse
# from django.http import HttpResponseRedirect
from django.http import JsonResponse
# from django.shortcuts import get_object_or_404
from django.shortcuts import render
# from django.urls import reverse
from django.utils import timezone
from django.views import generic

from .models import Protein, Result, Search, UserSearch
from .forms import ProtFindForm

# import json
# import simplejson
import bisect
import collections
import math
from pprint import pprint
import time
from typing import Deque, Dict, List, Set, Tuple


def make_pairs(flat_list: List) -> List[Tuple]:
    ''' return list of 2-tuples (pairs) from flat_list '''
    return list(zip(flat_list[0::2], flat_list[1::2]))


def make_trios(flat_list: List) -> List[Tuple]:
    ''' return list of 3-tuples (trios) from flat_list '''
    return list(zip(flat_list[0::3], flat_list[1::3], flat_list[2::3]))


def make_quads(flat: List) -> List[Tuple]:
    ''' return list of 4-tuples (quads) from flat_list '''
    return list(zip(flat[0::4], flat[1::4], flat[2::4], flat[3::4]))


def user_history(user: User, max_hist: int) -> List[Tuple]:
    '''
    Historical search results for one user.
    - Searches are in reverse chronological order (because ID increments).
    - Results parsed from CSV are in original found/saved order.
    - Uses make_pairs to parse the CSV, decoupled from DB string/CTE handling.
    '''
    return [(s.term, s.size, s.fmax, make_pairs(s.rcsv.split(',')))
            for s in Search.objects.filter(user=user)
                                   .order_by('-id')[:max_hist]]


def updated_or_created_term_result(acgt_str: str,
                                   max_show: int = 0) -> List[Tuple]:
    '''
    Update or create the saved Result for one search term.
    Return the list of (name, offset)-pairs.

    - Results form multiple searches and all users are combined.
    - Results parsed from CSV are in Protein-name-sorted order.
    - max_show parameter is ignored for now; 0 means no limit.

    We want a Django ORM equivalent of something like this SQL:

    SELECT DISTINCT INSTR(term, acgt_str) AS indx, term, GROUP_CONCAT(rcsv,':')
    FROM pfind_search
    WHERE indx > 0
    GROUP BY indx
    ORDER BY term;

    - Where acgt_str is the search-string argument to this function,
        and 'term' selects the field value to be searched.
    - We don't need or want to keep the 'term' field; only the
        substring offset 'indx' and the 'rcsv' value.
    - Using the DB's GROUP_CONCAT (or equivalent) function *might*
        be more efficient than processing multiple rows in Python,
        but we would need to split the concatenated results anyway,
        so let's not bother with an ORM-equivalent, and definitely
        don't invoke DB-specific functionality in a Func or Raw SQL
        statement.
    - Since we're not making the DB group results, we also don't
        gain anything by making the DB order them.  We need to
        process every row anyway.
    - Distinct, however, may be implemented more efficiently by the DB.

    Complexity:
    -   The operations of splitting, searching, joining,
        and rewriting search results are expected to be sublinear
        in time and linear in space on kilobytes,
    -   Substring searching on the Proteins is expected to be
        slightly sublinear or linear on megabytes.
    Thus it makes sense to update the results cache synchronously,
        whereas searching is dones asynchronously.
    '''
    qset = (Search.objects.annotate(indx=StrIndex('term', Value(acgt_str)))
                          .filter(indx__gt=0)
                          .values('indx', 'rcsv')
                          .distinct())
    pidmap = {}
    hitmap = collections.defaultdict(lambda: math.inf)
    for qres in qset:
        rcsv = qres['rcsv']
        if rcsv:
            toks = rcsv.split(",")
            indx = int(qres['indx'])
            for ncbi_name, found_idx in zip(toks[0::2], toks[1::2]):
                hitmap[ncbi_name] = min(hitmap[ncbi_name], int(found_idx) + indx)

    # List of the usual pairs (PName, MinIndex):
    pairs = [(name, hitmap[name]) for name in sorted(hitmap.keys())]

    # save results to cache
    term_rcsv = ','.join(','.join([name, str(indx)]) for name, indx in pairs)
    defaults = {'rcsv': term_rcsv, 'size': len(pidmap)}
    entry, created = Result.objects.update_or_create(term=acgt_str,
                                                     defaults=defaults)
    # server log:
    print('Cached Result for "%s": %s'
          % (acgt_str, ("CREATED" if created else "UPDATED")))

    return pairs


def update_saved_result(term_result: Result, user_pairs: List[Tuple]) -> None:
    '''
    Update an already saved Result row for one search term.

    Check results from one user search
    against those already saved from all users.
    If any result was not saved already, or has a smaller found offest,
    then update the saved Result.
    '''
    acgt_str = term_result.term
    term_rcsv = term_result.rcsv
    term_pairs = make_pairs(term_rcsv.split(','))   # sorted by name already
    term_names = [pair[0] for pair in term_pairs]
    do_updates = False
    name_index = 0
    user_pairs.sort(key=lambda x: x[0])             # sort by name
    for ncbi_name, found_idx in user_pairs:
        # Binary search saved Result for this found Protein name:
        name_index = bisect.bisect_left(term_names, ncbi_name, name_index)
        if (name_index < len(term_names)
                and term_names[name_index] == ncbi_name):
            # If ncbi_name was already found and saved in this Result, and
            # If the saved offset is greater than the newly found one,
            # replace it and flag the Result for updating:
            if (int(term_pairs[name_index][1]) > found_idx):
                found_idx_str = str(found_idx)
                term_pairs[name_index] = (ncbi_name, found_idx_str)
                do_updates = True
        else:
            # if the saved Result does not contain this Protein,
            # insert it at this index and flag the Result for updating:
            found_idx_str = str(found_idx)
            term_pairs.insert(name_index, (ncbi_name, found_idx_str))
            term_names.insert(name_index, ncbi_name)
            do_updates = True

    if do_updates:
        # NOTE: Docs say (get_or_create, update) is safer than (get, save).
        new_term_rcsv = ','.join(','.join(pair) for pair in term_pairs)
        Result.objects.filter(pk=acgt_str).update(date=timezone.now(),
                                                  size=len(term_pairs),
                                                  rcsv=new_term_rcsv)


@sync_to_async
def json_find_and_update(request) -> Dict:
    ''' AJAX wrapper to call dict_find_and_update '''
    # Make new search query.

    print("json_find_and_update:  request.GET: [{}]".format(request.GET))

    username = request.GET.get('username', '')
    acgt_str = request.GET.get('acgt_str', '')
    max_find = int(request.GET.get('max_find', 1))
    name_ord = request.GET.get('name_ord', False)
    max_hist = request.GET.get('max_hist', 0)

    print("json_find_and_update: {un:%s  ss:%s  mf:%s  no:%s  mh:%s}"
          % (username, acgt_str, max_find, name_ord, max_hist))

    user = User.objects.get(username=username)

    # # artificial delay for developing/debugging async:
    # secs = time.time()
    # print("json_find_and_update: %s at %d" % (username, secs))
    # wait = 12
    # time.sleep(wait)
    # secs = time.time()
    # print("json_find_and_update: %s at %d after %d"
    #       % (username, secs, wait))

    res_dict = dict_find_and_update(user,
                                    acgt_str,
                                    max_find,
                                    name_ord,
                                    max_hist)

    # NOTE: serializers.serialize('json', res_dict...) could also work.
    json_data = {
        'found_ct': res_dict['found_ct'],
        'dur_secs': res_dict['dur_secs'],
        'find_res': [(r['ncbi_name'], str(r['found_idx']))
                     for r in res_dict['find_res']]
    }
    response = JsonResponse(json_data, safe=True)
    return response


def dict_find_and_update(user: User,
                         acgt_str: str,
                         max_find: int,
                         name_ord: bool,
                         max_hist: int) -> Dict:
    '''
    - Search for proteins matching the search string
    - Update the Search and UserSearch objects
    - Return Python dictionary.
    Usage:
        - for synchronous Python, return the dict directly.
        - for AJAX, return dict wrapped in a JsonResponse
            - @see json_find_and_update

    The Raw Query is meant to be efficient, but is less portable.
    Using the Model objects / ORM might be just as good or better?

    Note from https://docs.djangoproject.com/:
        order_by('?') queries may be expensive and slow, depending on the
        database backend you’re using.
    But raw "SQL ORDER BY RANDOM()" is no better, because it probably *is*
        the backend translation.
    '''
    # Make new search query.
    res_dict = {}
    MAX_LIMIT = 1048576   # 2**20 is "enough for now"
    limit = max_find if max_find > 0 else MAX_LIMIT

    beg_time = time.perf_counter()
    # Note that StrIndex, like Sqlite's INSTR, uses 1-based indexing.
    ordered = 'pk' if name_ord else '?'
    find_res = (Protein.objects.annotate(
        found_idx=StrIndex('acgt_text', Value(acgt_str)) - 1)
        .filter(found_idx__gte=0)
        .order_by(ordered)
        .values('ncbi_name', 'found_idx')[:limit])
    dur_time = time.perf_counter() - beg_time
    dur_secs = "%6.2gs" % dur_time       # "%5.2g seconds"

    # Save Search and UserSearch rows in the DB
    found_ct = len(find_res)
    user_prs = []
    str_list = []
    for result in find_res:
        ncbi_name = result['ncbi_name']
        found_idx = result['found_idx']
        user_prs.append((ncbi_name, found_idx))
        str_list += [ncbi_name, str(found_idx)]
    user_rcsv = ','.join(str_list)
    search = Search.objects.create(user=user,
                                   size=found_ct,
                                   fmax=max_find,
                                   term=acgt_str,
                                   rcsv=user_rcsv)
    UserSearch.objects.create(user=user, search=search)

    # Get saved combined result for this term, if it exists;
    # else create it with this one new result -- which may be
    # empty if no Protein matched the search term.
    defaults = {'rcsv': user_rcsv}
    term_result, created = Result.objects.get_or_create(pk=acgt_str,
                                                        defaults=defaults)
    # If this Result was created here and now, we're done with it.
    # If it already existed, update it as needed.
    if not created and user_rcsv:
        update_saved_result(term_result, user_prs)

    res_dict['found_ct'] = found_ct
    res_dict['dur_secs'] = dur_secs
    res_dict['find_res'] = find_res
    return res_dict


def pfind_form(request):
    # initial dictionary overrides the form defaults:
    acgt_str = ""
    max_find = 1
    max_hist = 24       # TODO: parametrize?
    name_ord = False
    username = session_username = request.session.get('session_username', '')
    init_dict = {
        "acgt_str": acgt_str,
        "max_find": max_find,
        "name_ord": name_ord,
        "email_id": username,
    }
    form = ProtFindForm(request.POST or None, initial=init_dict)
    uhistory = []
    term_res = []
    user = User.objects.none()
    username = session_username = request.session.get('session_username', '')
    if username:
        try:
            user = User.objects.get(username=username)
            uhistory = user_history(user, max_hist)
        except Exception as err:
            print(err)

    # On submission of new form data:
    # - replace user, history
    # - execute (new) query
    if form.is_valid():
        acgt_str = form.cleaned_data.get("acgt_str")
        max_find = form.cleaned_data.get("max_find")
        name_ord = form.cleaned_data.get("name_ord")
        email_id = form.cleaned_data.get("email_id")

        # Get or create user (but don't use get_or_create):
        try:
            user = User.objects.get(username=email_id)
        except Exception as err:
            print(err)
            user = User.objects.create_user(
                username=email_id,
                email=email_id,
                password="")
            user.save()
        except IntegrityError as err:
            print(err)
            # TODO @sprax: display error message in form?  (Do not redirect)

        username = user.username
        request.session['session_username'] = username

        # If the user changed, get their history before executing new query.
        if username != session_username:
            uhistory = user_history(user, max_hist)

        # Find results saved from all user searches for this term.
        term_res = updated_or_created_term_result(acgt_str, max_hist)

        request.session.set_expiry(604800 * 4)  # remember sesskon for 4 weeks

    return render(request, 'pfind/pfind_form.html',
                  {'form': form,
                   'acgt_str': acgt_str,
                   'max_find': max_find,
                   'name_ord': name_ord,
                   'username': username,
                   'hist_len': len(uhistory),
                   'max_hist': max_hist,
                   'uhistory': uhistory,
                   'term_res': term_res,
                   'tres_len': len(term_res)})


class IndexView(generic.ListView):
    """ Use generic views: Less code is better """

    template_name = "pfind/list_view.html"
    context_object_name = "protein_list"

    def get_queryset(self, limit: int = 0):
        ''' Return all (or the last N) proteins, reverse-ordered by PK. '''
        if limit > 0:
            return Protein.objects.order_by("-pk")[:limit]
        return Protein.objects.order_by("-pk")

#
# def index(request):
#     p_list = Protein.objects.order_by('-pk')[:5]
#     output = ', '.join([p.ncbi_name for p in p_list])
#     return HttpResponse("This is the pfind index. <P> Some Protein names: "
#                         + output)


def detail(request, protein_id: int):
    try:
        protein = Protein.objects.get(pk=protein_id)
    except Protein.DoesNotExist:
        raise Http404("Protein %d not found" % protein_id)
    return render(request, "pfind/detail.html", {"protein": protein})


def ncbi_name(request, ncbi_name: str):
    try:
        protein = Protein.objects.get(ncbi_name=ncbi_name)
    except Protein.DoesNotExist:
        # NOTE: We could handle this with another template
        # instead of the 404, and set DEBUG = False in Django settings.
        raise Http404("Protein {} not found".format(ncbi_name))
    return render(request, "pfind/ncbi_name.html", {"protein": protein})
