# Generated by Django 3.1.3 on 2020-12-01 17:57

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pfind', '0004_search_size'),
    ]

    operations = [
        migrations.AddField(
            model_name='search',
            name='fmax',
            field=models.PositiveIntegerField(default=0),
        ),
    ]