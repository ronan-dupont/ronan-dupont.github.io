# Visitor counter

Activated on 14 September 2026. The footer displays a shared, estimated unique-visitor count across all pages and languages. It starts at activation; it does not reconstruct historical traffic. It estimates browser/network identities, not verified individual people.

## Implementation

`assets/js/visitor-counter.js` uses the public Busuanzi API at `https://bsz.iirose.cn/api`. The service stores the total and deduplicates visitors server-side. This small local client avoids loading third-party executable JavaScript. `_includes/visitor-counter.html` provides the EN/FR/JA label; `_sass/_visitor-counter.scss` styles it.

Only the canonical root `https://ronan-dupont.github.io/` is sent as the site identifier. The client sends no page path, query, title, referral URL or cookies. Like any HTTP server, the provider receives the request IP and browser User-Agent; its source hashes those for deduplication. A signed visitor marker is kept in local storage. A session marker avoids counting each navigation again.

Do Not Track and Global Privacy Control are respected. Requests run only on the production hostname. Timeout, blocked storage, unavailable service and malformed values are handled without showing fabricated numbers. The counter remains hidden when counting is unavailable or disabled.

## Service and accuracy

The free public host provides no availability or data-integrity guarantee. Shared networks, multiple devices, deleted storage and blockers can affect the estimate. For durable analytics, replace the provider with an owned analytics account or a self-hosted service.

## Sources

- [Public service and availability notice](https://busuanzi.9420.ltd/)
- [Server identity implementation](https://github.com/soxft/busuanzi/blob/main/app/middleware/identity.go)
- [Shared count implementation](https://github.com/soxft/busuanzi/blob/main/core/count.go)
- [GET and POST behavior](https://github.com/soxft/busuanzi/blob/main/app/controller/api.go)

The local client was checked for production scoping, shared-count rendering, root-only requests, session reads, privacy preferences, malformed data and network failures. CORS and read-only API access were verified against the real service.
