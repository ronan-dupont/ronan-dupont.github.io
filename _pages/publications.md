---
layout: archive
title: "Publications"
permalink: /publications/
lang: en
translation_key: publications
author_profile: true
---

{% if site.author.googlescholar %}
  You can also find my articles on <u><a href="{{site.author.googlescholar}}">my Google Scholar profile</a>.</u>
{% endif %}

{% include base_path %}

<p class="archive-intro">Peer-reviewed articles, conference papers, manuscripts and selected academic reports.</p>

{% for post in site.publications reversed %}
  {% include archive-single-localized.html %}
{% endfor %}
