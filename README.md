# Protein Finder

Finds **exact** matches of DNA sequences (such as 'acgttgca') in a database of proteins.

## Form Inputs:
DNA string, max number of matches, user email, and find/display order (name order v. random)

## Output Tables:
- New Results: matches found from just-submitted form.  
    - A match comprises an NCBI Protein name and the zero-base offset of the first known occurence of the search term.
- User Search History: previous searches by the same user, shown in reverse chronological order.
- Saved Term Results: previous search results for the current search term, aggregated from all user's searches, always shown in NCIB-name order.
    - Saved results include substring matches: If 'acgtacgt' was found at offset X, then the substring 'tacg' was found at offset X + 3.
    - Any offset shown is the smallest known.  A later search for the substring may find a smaller offset.

Each new search is asynchronous, so the "New" results dispplayed are those most recently returned, not necessary those of the most recently submitted search.  If one search is leap-frogged by another, its results show up in the Search History and Saved Resulta tables rather than in New Results.

Users are identified by email only, which is stored in a cookie for session continuity.  There are no User passwords.  The Admin password is an 8-letter local name for I-90 in PascalCase.

## Protein Data Source

DNA sequence data was download from NCBI Entrez using various Python scripts and scraping methods.
The simplest way was to use the eponymous script from https://github.com/kblin/ncbi-acc-download.
Example:

```shell
ncbi-acc-download NC_027867 --format fasta -v
```
Then it is easy to use a shell command such as sed, or an editor such as vim,
to convert the downloaed text into a single string containing only lowercase
letters 'a', 'c', 'g', 't', or in some cases 'n' (to wit, from NC_008724 and
NC_023640).

## Notes

Note that partial alignments, such as found using:

    from Bio.Seq import Seq
    alignments = Bio.pairwise2.align.globalxx(Seq("ACGTTGCA"), protein_seq)

are not supported by this site.

### Punch List:

- ~Cookies to resume session~
- Done: ~Async find and page updates~
    - Not found: ~Use newer Django built-in or contrib~
    - Done: ~Use AJAX/JQuery if this newer stuff fails right out of the box, or, as it seems now, is not meant to replace AJAX all at once.~
    - ~Test using artificial delay.~
- Nice To Have: Done: ~Aggregate saved results from all users.~
- Nice To Have: Not warrented now: ~In-memory cache with warming and persistence.~

### Dev-ony NOTES (down in the weeds):

```python
>>> len([(s.id, s.term, s.date.minute, s.rcsv.split(',')) for s in Search.objects.filter(term__contains='gtca').order_by('user').distinct().reverse()])
8
>>> len([(s.id, s.term, s.date.minute, s.rcsv.split(',')) for s in Search.objects.filter(term__contains='ggtca').order_by('user').distinct().reverse()])
4
```

The 8-letter admin password is not ezpass, but where it's used.

```python
>>> len([(s.id, s.term, s.date.minute, s.rcsv.split(',')) for s in Search.objects.filter(term__contains='gtca').order_by('user').distinct().reverse()])
8
>>> len([(s.id, s.term, s.date.minute, s.rcsv.split(',')) for s in Search.objects.filter(term__contains='ggtca').order_by('user').distinct().reverse()])
4
```

More nice-to-have's: [Nice To Have] Cache previous searches in-memory and in DB.
- On server startup, load the cache from the DB.
- If the search string is a key in the in-memory cache,
    retrieive that.
    Otherwise do the DB raw query and cache the results.
- The DB relations are:
    - User -> many Querys (history)
    - Query -> many Results (name, index)-pairs
- Hit the Protein table only when a Query is not
    in the cache?
- Invalidate (parts of) the cache to force new searches.
- Optimize with a link table?
    - Yes, not for performance on a small scale,
        but for cleaner design/sepration of concerns.
        So: [ UserId | Query ]

### More on saving/caching results from one or all users

#### Two main options -- parse and rewrite CTE results inside or outside the DB:

1.  Using a View or Table: create the view with a selection such as

    SELECT term, rcsv FROM pfind_search
    WHERE  term LIKE ('%'||'acgt'||'%')

for every distinct search term (replacing the literal 'acgt' above),
regardless of user.  But we need the offset index of a search term
if it is a substring of a larger search term, so we may need to do
a correlated subquery instead of using EXISTS on a simple subquery,
or create a (temporary) view for a WITH-clause.

* If we wanted a user-specific cache, we'd either:
        - include a user ID in each Search row, or,
        do an inner join of User and Search on user ID
        via UserSearch.

SIMPLEST way, no tracking and inference sub-search-strings:
sqlite> SELECT term, GROUP_CONCAT(rcsv,':') FROM pfind_search GROUP BY
        term ORDER BY term;
(And do the rest in Python, avoiding CSV parsing in SQLite.
PostgreSQL would be more plausible, but still limits flexibility)

sqlite> SELECT id, date, term, INSTR(term, 'acgt') - 1, rcsv
        FROM pfind_search
        WHERE term LIKE '%acgt%'
        ORDER BY rcsv;

sqlite> SELECT term, GROUP_CONCAT(rcsv,',') FROM pfind_search WHERE
        EXISTS (SELECT pfs.term, instr(pfs.term, term) AS idx FROM
        pfind_search AS pfs WHERE pfs.term LIKE ('%'||term||'%'))
        GROUP BY term;

sqlite> SELECT term, GROUP_CONCAT(rcsv,',') FROM pfind_search WHERE
        EXISTS (SELECT pfs.term, instr(pfs.term, term) AS off FROM
        pfind_search AS pfs WHERE off > 0) GROUP BY term;

#### Long, concatenated results:
sqlite> select distinct term, GROUP_CONCAT(rcsv, ':') from pfind_search
        where exists (select distinct pfs.term from pfind_search AS pfs
        where instr(pfs.term, term) > 0) group by term limit 7;


#### One Way (without creating permanent view of unique search terms):
Use WITH-clause to create temporary table U of unique search terms
from the whole Search table pfind_search as S;
and for each U.term, find all previous matches of all terms of which
U.term is a substring, along with the index of U.term in S.term.
Example: Any Protein sequence match for 'acgtacgt' is also a match
for 'gtac', a substring of 'acgtacgt' at index 2.  So if 'acgtacgt'
was found in P at offset X, then in effecct, 'gtac' was also found,
but at offset X + 2.  Of course 'gtac' may also appear earlier in P,
and we would have already found it there if we had ever explicitly
searched for 'gtac' in P.  But we may have searched for 'gtac' in
some limited capacity (as in max_find != 0), and so 'gtac' exists
as a previus search term without its having been applied to P.
So we take the minimum offset of 'gtac' in P that we have already
found.  This offset may be preempted by a new search.

2. More portable: Use mostly Python to read and re-write cache...

DONE @sprax: Re-order the saved Result table operations:
- Update the Result table with the results of each new search;
    thus only at the end of *find_and_update.
- On valid form submission, parse and render the contents
    already stored in the Result table.
- [Nice To Have] It could be even more responsive to update the
    Result table in a separate thread spun off by *find_and_update,
    but the Django ORM may not be sufficiently thread-safe.
