---
layout: archive
title: Contact
permalink: /contact/
lang: en
translation_key: contact
author_profile: false
profile_page: true
---

<div class="profile-page profile-page--contact">
  {% include contact-identity.html eyebrow="Get in touch" invitation="Let’s talk about research." description="You can contact me about applied mathematics, scientific computing and teaching." portrait_alt="Illustrated portrait of Ronan Dupont" %}
  <div class="profile-page__contact-grid">
    <div class="profile-page__contact-methods">
      <a class="profile-page__contact-link" href="mailto:r-dupont@na.nuap.nagoya-u.ac.jp">
        <span class="profile-page__contact-icon"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><rect x="3" y="5" width="18" height="14" rx="3"/><path d="m4 7 8 6 8-6"/></svg></span>
        <span><span class="profile-page__contact-label">Email</span><strong>r-dupont@na.nuap.nagoya-u.ac.jp</strong></span>
        <span class="profile-page__contact-arrow" aria-hidden="true">↗</span>
      </a>
      <a class="profile-page__contact-link" href="tel:+33670908808">
        <span class="profile-page__contact-icon"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><path d="M8 3H5a2 2 0 0 0-2 2c0 8.8 7.2 16 16 16a2 2 0 0 0 2-2v-3l-5-2-2 2a13 13 0 0 1-6-6l2-2-2-5Z"/></svg></span>
        <span><span class="profile-page__contact-label">Telephone</span><strong>+33 6 70 90 88 08</strong></span>
        <span class="profile-page__contact-arrow" aria-hidden="true">↗</span>
      </a>
    </div>
    <section class="profile-page__location" aria-labelledby="contact-location-heading">
      <p class="profile-page__eyebrow"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" aria-hidden="true"><path d="M19 10c0 5-7 11-7 11S5 15 5 10a7 7 0 1 1 14 0Z"/><circle cx="12" cy="10" r="2.5"/></svg>Based in Nagoya, Japan</p>
      <h2 id="contact-location-heading">Nagoya University</h2>
      <p>Zhang–Sogabe Laboratory<br>Graduate School of Engineering</p>
      <p class="profile-page__location-role">Postdoctoral researcher in numerical linear algebra</p>
      <a class="profile-page__text-link" href="{{ '/cv/' | relative_url }}">View my academic profile <span aria-hidden="true">→</span></a>
    </section>
  </div>
  <section class="profile-page__profiles" aria-labelledby="contact-profiles-heading">
    <div class="profile-page__section-heading"><h2 id="contact-profiles-heading">Research and professional profiles</h2><p>Find my researcher identifiers, publications and code through these profiles.</p></div>
    <div class="profile-page__profile-links">
      <a href="{{ site.author.orcid | escape }}"><span><strong>ORCID</strong><span>Researcher identifier</span></span><span aria-hidden="true">↗</span></a>
      <a href="{{ site.author.researchgate | escape }}"><span><strong>ResearchGate</strong><span>Research profile</span></span><span aria-hidden="true">↗</span></a>
      <a href="https://github.com/{{ site.author.github | escape }}"><span><strong>GitHub</strong><span>Code and repositories</span></span><span aria-hidden="true">↗</span></a>
      <a href="https://www.linkedin.com/in/{{ site.author.linkedin | escape }}"><span><strong>LinkedIn</strong><span>Professional profile</span></span><span aria-hidden="true">↗</span></a>
    </div>
  </section>
  <nav class="profile-page__languages profile-page__languages--footer" aria-label="Contact page languages">
    <a href="{{ '/contact/' | relative_url }}" lang="en" hreflang="en"{% if page.lang == 'en' %} aria-current="page"{% endif %}>English</a>
    <a href="{{ '/fr/contact/' | relative_url }}" lang="fr" hreflang="fr"{% if page.lang == 'fr' %} aria-current="page"{% endif %}>Français</a>
    <a href="{{ '/ja/contact/' | relative_url }}" lang="ja" hreflang="ja"{% if page.lang == 'ja' %} aria-current="page"{% endif %}>日本語</a>
  </nav>
</div>
