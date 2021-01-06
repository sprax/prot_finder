# @file: pfind/urls.py
'''
NOTE: To handle overloaded URLs properly,
the <int:protein_id> path must precede the <str:ncbi_name> path in the list.
'''
from django.conf.urls import url
from django.contrib import admin
from django.urls import path, include
from . import views

app_name = "pfind"

urlpatterns = [
    path("", views.pfind_form, name="index"),
    # path("find/", views.pfind_form, name="pfind_form"),
    # path("list/", views.IndexView.as_view(), name="list"),
    path("<str:ncbi_name>/", views.ncbi_name, name="ncbi_name"),
    url(r'^ajax/json_find_and_update/$',
        views.json_find_and_update, name='json_find_and_update'),
]
