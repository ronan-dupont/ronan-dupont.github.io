---
layout: archive
title: Contact
permalink: /fr/contact/
lang: fr
translation_key: contact
author_profile: false
profile_page: true
---

<div class="profile-page profile-page--contact">
  <section class="profile-page__contact-hero" aria-labelledby="contact-heading">
    <p class="profile-page__eyebrow">Échanger</p>
    <h2 id="contact-heading" class="profile-page__lead-title">Parlons de recherche.</h2>
    <p class="profile-page__lead">Vous pouvez me contacter au sujet des mathématiques appliquées, du calcul scientifique ou de l’enseignement.</p>
  </section>
  <div class="profile-page__contact-grid">
    <div class="profile-page__contact-methods">
      <a class="profile-page__contact-link" href="mailto:r-dupont@na.nuap.nagoya-u.ac.jp">
        <span class="profile-page__contact-icon"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><rect x="3" y="5" width="18" height="14" rx="3"/><path d="m4 7 8 6 8-6"/></svg></span>
        <span><span class="profile-page__contact-label">E-mail</span><strong>r-dupont@na.nuap.nagoya-u.ac.jp</strong></span>
        <span class="profile-page__contact-arrow" aria-hidden="true">↗</span>
      </a>
      <a class="profile-page__contact-link" href="tel:+33670908808">
        <span class="profile-page__contact-icon"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><path d="M8 3H5a2 2 0 0 0-2 2c0 8.8 7.2 16 16 16a2 2 0 0 0 2-2v-3l-5-2-2 2a13 13 0 0 1-6-6l2-2-2-5Z"/></svg></span>
        <span><span class="profile-page__contact-label">Téléphone</span><strong>+33 6 70 90 88 08</strong></span>
        <span class="profile-page__contact-arrow" aria-hidden="true">↗</span>
      </a>
    </div>
    <section class="profile-page__location" aria-labelledby="contact-location-heading">
      <p class="profile-page__eyebrow"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" aria-hidden="true"><path d="M19 10c0 5-7 11-7 11S5 15 5 10a7 7 0 1 1 14 0Z"/><circle cx="12" cy="10" r="2.5"/></svg>À Nagoya, au Japon</p>
      <h2 id="contact-location-heading">Université de Nagoya</h2>
      <p>Laboratoire Zhang–Sogabe<br>Graduate School of Engineering</p>
      <p class="profile-page__location-role">Chercheur postdoctoral en algèbre linéaire numérique</p>
      <a class="profile-page__text-link" href="{{ '/cv-fr/' | relative_url }}">Consulter mon parcours <span aria-hidden="true">→</span></a>
    </section>
  </div>
  <section class="profile-page__profiles" aria-labelledby="contact-profiles-heading">
    <div class="profile-page__section-heading"><h2 id="contact-profiles-heading">Profils scientifiques et professionnels</h2><p>Retrouvez mes identifiants de chercheur, mes publications et mon code sur ces profils.</p></div>
    <div class="profile-page__profile-links">
      <a href="{{ site.author.orcid | escape }}"><span><strong>ORCID</strong><span>Identifiant de chercheur</span></span><span aria-hidden="true">↗</span></a>
      <a href="{{ site.author.researchgate | escape }}"><span><strong>ResearchGate</strong><span>Profil scientifique</span></span><span aria-hidden="true">↗</span></a>
      <a href="https://github.com/{{ site.author.github | escape }}"><span><strong>GitHub</strong><span>Code et dépôts</span></span><span aria-hidden="true">↗</span></a>
      <a href="https://www.linkedin.com/in/{{ site.author.linkedin | escape }}"><span><strong>LinkedIn</strong><span>Profil professionnel</span></span><span aria-hidden="true">↗</span></a>
    </div>
  </section>
  <nav class="profile-page__languages profile-page__languages--footer" aria-label="Langues de la page contact">
    <a href="{{ '/contact/' | relative_url }}" lang="en" hreflang="en"{% if page.lang == 'en' %} aria-current="page"{% endif %}>English</a>
    <a href="{{ '/fr/contact/' | relative_url }}" lang="fr" hreflang="fr"{% if page.lang == 'fr' %} aria-current="page"{% endif %}>Français</a>
    <a href="{{ '/ja/contact/' | relative_url }}" lang="ja" hreflang="ja"{% if page.lang == 'ja' %} aria-current="page"{% endif %}>日本語</a>
  </nav>
</div>
