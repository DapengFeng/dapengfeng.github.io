// Use Playwright's pinned Chromium locally and in CI. An alternate browser must
// be selected explicitly so an installed desktop browser cannot mask CI failures.
export function chromiumExecutable() {
 return process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE || undefined;
}
