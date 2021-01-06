'''
@file: models.py

SCHEMA
Protein     proteins as DNA sequences (PK = NCBI_name, ACGT-sequence)
User        users who search proteins [standard Django]
Search      particular searches (id, date, count, term, result-as-CTE)
UserSearch  each row links one User ID to one Search ID
Result      all known results for each search term (PK=search_term, date, CSV)

Enables two search result views: Historical or Aggregated.
Historical: list of searchs by user, sorted by date or term.
            The same search term may appear many times with same or
            different results.
Aggregated: All known results for each search term,
            abstracted from all users and dates,
            but holding a last-updated date (and possibly
            the ID of the last User to search for the term).
'''
# from datetime import datetime
# from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone


class Protein(models.Model):
    '''
    Just the ID (implicit), name, and sequence.
    'pub_date' could of use if new proteins are added
    often, and search were so expensive that we might
    want to search only those added since the date of
    some previous search.

    TODO @sprax: Make ncbi_name the PK?
    '''
    # NOTE: names are assumed to be alphanumeric.  Underscore is find,
    #       not punctuation that would need to be escaped in CSV format.
    ncbi_name = models.CharField(max_length=10, primary_key=True)
    acgt_text = models.CharField(max_length=3200000)
    # pub_date = models.DateTimeField('date published')


class Search(models.Model):
    '''
    - Search term,
    - DateTime of search,
    - Result of single search, stored as an unsized CSV list
        of (name, index)-pairs, sorted by name.  May be treated as a CTE.

    Examples:
    term:   "acgtaccattag"
    rcsv:   "NCBI_nameA,offsetA,NCBI_nameB,offsetB"
            where NCBI_name is a 9-character ID such as NC_123456
            and offsets are to be converted to non-negative integers.
    - We could remove the 'user' ForeignKey and rely only
        on the UserSearch link table.
    '''
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    date = models.DateTimeField(default=timezone.now)
    size = models.PositiveIntegerField(default=0)   # num found matches
    fmax = models.PositiveIntegerField(default=0)   # max_find (limit)
    term = models.CharField(max_length=4092)        # saerch term
    rcsv = models.CharField(max_length=16368)       # comma-separated results


class UserSearch(models.Model):
    '''
    UserSearch  each row links one User ID to one Search ID.
    - This is somewhat redundant if Search includes a User ID.
    '''
    user = models.ForeignKey(User, null=False,
                             on_delete=models.CASCADE)
    search = models.ForeignKey(Search, null=False,
                               on_delete=models.CASCADE)
    date = models.DateTimeField(default=timezone.now)


class Result(models.Model):
    '''
    - Search term  (unique key),
    - DateTime of most recent search == last update,
    - Results accumulated from multiple searches,
        stored as an unsized, sorted CSV list, keeping only the minimum
        offset index of the term found (so far) for each Protein sequence.

    Example:
    term:   "acgtaccattag"
    rcsv:   "NCBI_nameA,offsetA,NCBI_nameB,offsetB"
            where NCBI_name is a 9-character ID such as NC_123456
            and offsets are to be converted to non-negative integers.
    '''
    date = models.DateTimeField(default=timezone.now)
    size = models.PositiveIntegerField(default=0)   # num found matches
    term = models.CharField(max_length=4092, primary_key=True)  # saerch term
    rcsv = models.CharField(max_length=16368)       # comma-separated results
