---
layout: page
permalink: /publications/
title: Publications
description: Scientific publications in international journals and conferences.
years_papers: [2017,2020]
years_papers: [2016,2019]
nav: true
---

Work in progress

### International Peer-Reviewed Journals
<div class="publications">

{% for y in page.years %}
  <h2 class="year">{{y}}</h2>
  {% bibliography -f papers -q @*[year={{y}}]* %}
{% endfor %}

</div>

### Conferences
<div class="publications">

{% for y in page.years %}
  <h2 class="year">{{y}}</h2>
  {% bibliography -f conferences -q @*[year={{y}}]* %}
{% endfor %}

</div>
