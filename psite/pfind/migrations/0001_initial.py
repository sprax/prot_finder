# Generated by Django 3.1.3 on 2020-11-19 20:37

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('ncbi_name', models.CharField(max_length=10)),
                ('acgt_text', models.CharField(max_length=2000000)),
            ],
        ),
    ]