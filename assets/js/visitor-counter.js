/* Site-wide estimated visitors. The count is shared by the remote service.
 * Only the canonical site root is sent; no page URL, title or referrer is sent.
 * Local storage holds only the service's signed visitor marker for deduplication.
 * No third-party JavaScript is loaded and no cookies are sent. */
(function () {
  'use strict';
  var counter = document.querySelector('[data-visitor-counter]');
  if (!counter || window.location.hostname !== 'ronan-dupont.github.io') return;
  if (navigator.doNotTrack === '1' || window.doNotTrack === '1' || navigator.globalPrivacyControl) return;
  if (typeof window.fetch !== 'function' || typeof window.AbortController !== 'function') return;

  var endpoint = 'https://bsz.iirose.cn/api';
  var siteRoot = 'https://ronan-dupont.github.io/';
  var identityKey = 'rd-visitor-identity';
  var sessionKey = 'rd-visitor-counted';
  var headers = { 'x-bsz-referer': siteRoot };
  var counted = false;
  try {
    var identity = window.localStorage.getItem(identityKey);
    if (identity && /^[a-f0-9]{32}\.[a-f0-9]{40,64}$/i.test(identity)) headers.Authorization = 'Bearer ' + identity;
  } catch (_) { /* Private browsing may disable local storage. */ }
  try { counted = window.sessionStorage.getItem(sessionKey) === '1'; } catch (_) {}

  var abort = new AbortController();
  var timeout = window.setTimeout(function () { abort.abort(); }, 8000);
  window.fetch(endpoint, {
    method: counted ? 'GET' : 'POST',
    headers: headers,
    credentials: 'omit',
    referrerPolicy: 'no-referrer',
    cache: 'no-store',
    signal: abort.signal
  }).then(function (response) {
    if (!response.ok) throw new Error('Counter unavailable');
    var identity = response.headers.get('Set-Bsz-Identity');
    if (identity && /^[a-f0-9]{32}\.[a-f0-9]{40,64}$/i.test(identity)) {
      try { window.localStorage.setItem(identityKey, identity); } catch (_) {}
    }
    return response.json();
  }).then(function (result) {
    var count = result && result.success && result.data && result.data.site_uv;
    if (typeof count !== 'number' || !Number.isSafeInteger(count) || count < 0) return;
    counter.querySelector('[data-visitor-count]').textContent = count.toLocaleString(document.documentElement.lang || 'en');
    counter.hidden = false;
    try { window.sessionStorage.setItem(sessionKey, '1'); } catch (_) {}
  }).catch(function () {
    // Leave the small counter hidden when blocked, offline or unavailable.
  }).then(function () { window.clearTimeout(timeout); });
})();
