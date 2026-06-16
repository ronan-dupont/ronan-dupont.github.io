---
layout: archive
title: "Publications"
permalink: /fr/publications/
lang: fr
translation_key: publications
author_profile: true
---

{% if site.author.googlescholar %}
  Vous pouvez également retrouver mes articles sur <u><a href="{{site.author.googlescholar}}">mon profil Google Scholar</a>.</u>
{% endif %}

{% include base_path %}

<p class="archive-intro">Articles évalués par les pairs, actes de conférence, manuscrits et rapports académiques sélectionnés.</p>

{% for post in site.publications reversed %}
  {% include archive-single-localized.html %}
{% endfor %}
