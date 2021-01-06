

# import datetime
# from django.utils import timezone

from django.test import TestCase, RequestFactory, AsyncRequestFactory

from pfind.models import Protein

from django.shortcuts import render
from django.http import HttpRequest
from django.template.response import TemplateResponse

import collections
import math


class ProteinModelTests(TestCase):
    def test_protein_crud(self):
        """
        CRUD
        """
        ncbi_name = "NC_test"
        protein = Protein(ncbi_name=ncbi_name)
        self.assertIs(protein.ncbi_name, ncbi_name)

        acgt_text = "ACGTTGCA"
        protein.acgt_text = acgt_text
        self.assertEqual(protein.acgt_text, acgt_text)


def my_template_view(request, template_name):
    """
    Here's a view using TemplateResponse. This test will run
    without needing to set up urls or add templates to the file system.
    """
    return TemplateResponse(request, template_name, {"foo": "bar"})


class TemplateResponseTest(TestCase):
    def test_get(self):
        request = RequestFactory().get("/")
        response = my_template_view(request, "template.html")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context_data["foo"], "bar")


def reduce_result_cte_to_csv(oldcsv: str, offset: int) -> str:
    '''
    - Replace duplicate NCBI Protein names
        by the one instance with the minimum offset;
    - Return a new CSV string that is sorted by Protein name
    (not by the DB's ID/PK for the row in Protein.)
    '''
    pidmap = {}
    hitmap = collections.defaultdict(lambda: math.inf)
    result = oldcsv.split(":")
    for row in result:
        tks = row.split(",")
        for p_id, ncbi_name, found_idx in zip(tks[0::3], tks[1::3], tks[2::3]):
            pidmap[ncbi_name] = p_id
            hitmap[ncbi_name] = min(hitmap[ncbi_name], int(found_idx) + offset)
    fields = []
    for name in sorted(hitmap.keys()):
        fields.append(pidmap[name])
        fields.append(name)
        fields.append(str(hitmap[name]))
    return ','.join(fields)


class ResultModelTests(TestCase):
    def test_result_reduce(self):
        """
        Expect duplicate NCBI Protein names to be replaced
        by the one instance with the minimum offset, and that
        new CSV string will be sorted by Protein name
        (not by the DB's ID/PK for the row in Protein.)

        interm is a proper substring of exterm, starting at offset.
        oldcsv represents the overlapping search results from searching
        for exterm 3 times.
        Note that if exterm was found at offset X,
        then interm can be found at index X + offset,
        but this index need not be the global minimum.
        """

        interm = "cgtaa"
        offset = 2
        exterm = "tccgtaaat"
        self.assertLessEqual(offset + len(interm), len(exterm))
        self.assertEqual(interm, exterm[offset:offset + len(interm)])
        oldcsv = ("1,NC_000111,888888,2,NC_000002,23456,3,NC_008333,2222:"
                  "1,NC_000111,555555,3,NC_008333,27672,5,NC_054321,5678:"
                  "2,NC_000002,111111,4,NC_444444,12345,9,NC_987654,9988"
                  )
        name_0 = "NC_000002"
        name_9 = "NC_987654"
        expect = ("2,NC_000002,23458,1,NC_000111,555557,3,NC_008333,2224,"
                  "5,NC_054321,5680,4,NC_444444,12347,9,NC_987654,9990"
                  )

        # redundant code check for reduce_result_cte_to_csv
        pidmap = {}
        hitmap = collections.defaultdict(lambda: math.inf)
        result = oldcsv.split(":")
        self.assertEqual(len(result), 3)
        for row in result:
            tks = row.split(",")
            self.assertEqual(len(tks), 9)
            for p_id, ncbi_name, found_idx in zip(tks[0::3], tks[1::3], tks[2::3]):
                pidmap[ncbi_name] = p_id
                hitmap[ncbi_name] = min(hitmap[ncbi_name], int(found_idx) + offset)
        self.assertEqual(len(hitmap), 6)
        fields = []

        sorted_names = list(sorted(hitmap.keys()))
        self.assertEqual(len(sorted_names), 6)
        self.assertEqual(name_0, sorted_names[0])
        self.assertEqual(name_9, sorted_names[-1])

        for name in sorted_names:
            fields.append(pidmap[name])
            fields.append(name)
            fields.append(str(hitmap[name]))
        newcsv = ','.join(fields)
        print("\t newcsv", newcsv)
        print("\t expect", expect)
        self.assertEqual(newcsv, expect)

        reduced = reduce_result_cte_to_csv(oldcsv, offset)
        self.assertEqual(reduced, expect)
