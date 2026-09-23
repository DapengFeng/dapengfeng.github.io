import {existsSync} from 'node:fs';

// Honor an explicit executable; otherwise use the desktop browser when available,
// or let Playwright select the version installed by CI.
export function chromiumExecutable() {
 const explicit=process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE;
 if(explicit)return explicit;
 const desktop='/home/jarvis/.cache/ms-playwright/chromium-1234/chrome-linux64/chrome';
 return existsSync(desktop)?desktop:undefined;
}
