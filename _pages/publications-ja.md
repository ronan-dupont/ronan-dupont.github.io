---
layout: archive
title: "研究業績"
permalink: /ja/publications/
lang: ja
translation_key: publications
author_profile: true
---

{% if site.author.googlescholar %}
  論文は <u><a href="{{site.author.googlescholar}}">Google Scholar profile</a></u> でも確認できます。
{% endif %}

{% include base_path %}

<p class="archive-intro">査読付き論文、会議論文、学位論文、主要な学術レポートをまとめています。</p>

{% for post in site.publications reversed %}
  {% include archive-single-localized.html %}
{% endfor %}
