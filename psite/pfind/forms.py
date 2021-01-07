# @file: forms.py

from django import forms


class ProtFindForm(forms.Form):
    '''
    acgt_str    DNA sequence string to search for in model's acgt_text field.
    name_ord    Search proteins in roughly the sort order of
                their names, even though in parallel searches
                may make find order indeterminate.  Mutliple
                matches will be presented in name-sorted order.
    max_find    Limit the search to this number of matches.
                0 means no limit.
                1 is the default; quit at the first match.
                This is a full stop, not a pause,
                There is no "resume from last found match."
    email_id    Token to identify a returning user.
                Validated but not authenticated.
    '''
    acgt_str = forms.CharField(
        label="DNA_str",
        min_length=3,
        initial="acgt",     # override this to prevent accidental searche
        help_text="DNA sequence to find (length > 2); example: 'gattattaccat'")
    max_find = forms.IntegerField(
        label="Limited",
        initial=1,
        help_text="Limit the search to N matches; 0 means no limit.")
    email_id = forms.EmailField(
        label="EmailID",
        help_text="Email address used to identify a returning user.")
    name_ord = forms.BooleanField(
        label='Ordered',
        initial=False,
        required=False,
        help_text=("Find and show matches in name-sorted order "
                   " v. random find order."))

    def clean_acgt_str(self):
        acgt_str = self.cleaned_data.get("acgt_str")
        if any(c not in "acgt" for c in acgt_str):
            raise forms.ValidationError(
                "The 'DNA string' must contain only these letters: a c g t")
        return acgt_str

    def clean_max_find(self):
        min_find = 0
        max_find = self.cleaned_data.get("max_find")
        if max_find < min_find:
            raise forms.ValidationError(
                "The 'Max Find Limit' must be at least %d" % min_find)
        return max_find
