// Keep previews and CI out of production statistics, and respect privacy signals.
const config = document.querySelector('script[data-analytics-token][data-analytics-host]');
const optedOut = navigator.globalPrivacyControl === true || navigator.doNotTrack === '1' || window.doNotTrack === '1';
if (config && location.protocol === 'https:' && location.hostname === config.dataset.analyticsHost && !optedOut) {
  const beacon = document.createElement('script');
  beacon.type = 'module';
  beacon.src = 'https://static.cloudflareinsights.com/beacon.min.js';
  beacon.setAttribute('data-cf-beacon', JSON.stringify({token: config.dataset.analyticsToken}));
  document.body.append(beacon);
}
