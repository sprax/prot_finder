# Generated by Django 3.1.3 on 2020-12-24 06:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pfind', '0007_auto_20201221_1531'),
    ]

    operations = [
        migrations.AlterField(
            model_name='protein',
            name='ncbi_name',
            field=models.CharField(max_length=10, unique=True),
        ),
    ]