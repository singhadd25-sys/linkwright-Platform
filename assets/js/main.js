(() => {
  const modal = document.getElementById("auth-modal");
  const signupPanel = document.getElementById("signup-panel");
  const loginPanel = document.getElementById("login-panel");
  const tabs = [...document.querySelectorAll("[data-auth-tab]")];
  function showAuth(mode) {
    modal.classList.add("open");
    modal.setAttribute("aria-hidden", "false");
    setMode(mode);
  }
  function closeAuth() {
    modal.classList.remove("open");
    modal.setAttribute("aria-hidden", "true");
  }
  function setMode(mode) {
    const signup = mode !== "login";
    signupPanel.hidden = !signup;
    loginPanel.hidden = signup;
    tabs.forEach(tab => tab.classList.toggle("active", tab.dataset.authTab === (signup ? "signup" : "login")));
  }
  const API = ((window.LINKWRIGHT_CONFIG || {}).API_BASE_URL || "").replace(/\/+$/, "");
  // Accounts, OTP email and sessions live on the HTTPS backend (same-origin cookies), so the
  // public site sends visitors there instead of posting credentials cross-origin.
  document.querySelectorAll("[data-open-auth]").forEach(button => button.addEventListener("click", () => {
    if (API) { window.location.href = API + (button.dataset.openAuth === "login" ? "/" : "/signup"); return; }
    showAuth(button.dataset.openAuth);
  }));
  tabs.forEach(tab => tab.addEventListener("click", () => setMode(tab.dataset.authTab)));
  document.querySelector(".modal-close").addEventListener("click", closeAuth);
  modal.addEventListener("click", event => { if (event.target === modal) closeAuth(); });
  document.addEventListener("keydown", event => { if (event.key === "Escape") closeAuth(); });
  document.getElementById("year").textContent = new Date().getFullYear();

  if (!API) {
    document.getElementById("signup-status").textContent = "Preview only: set API_BASE_URL in assets/js/config.js to the deployed Linkwright backend to enable registration.";
    document.getElementById("login-status").textContent = "Preview only: no credentials are sent or stored by this page.";
    ["signup-form","login-form"].forEach(id => document.getElementById(id).addEventListener("submit", e => e.preventDefault()));
  }
  const mobile = document.querySelector(".mobile-menu");
  mobile.addEventListener("click", () => {
    const nav = document.querySelector(".desktop-nav");
    const open = nav.style.display === "flex";
    nav.style.display = open ? "none" : "flex";
    if (!open) { nav.style.position="absolute"; nav.style.top="72px"; nav.style.left="0"; nav.style.right="0"; nav.style.padding="20px"; nav.style.background="#0c1624"; nav.style.flexDirection="column"; nav.style.border="1px solid #26384d"; }
  });
})();