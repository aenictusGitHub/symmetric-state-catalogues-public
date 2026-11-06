/* Typeset state picker. The original select remains the data/change interface
 * and is the fallback when JavaScript is unavailable. */
(function () {
  'use strict';

  let select, wrapper, trigger, display, menu, currentState = null;
  let emptyText = 'Explore a point on the plot', active = -1, options = [], optionNodes = [];
  let optionSignature = '', search = '', searchTime = 0;
  const originalLabels = new WeakMap();
  const plainNames = {
    'n4-tetra': 'Tetrahedral state',
    'n5-psia': 'psi 5 (a)', 'n5-psib': 'psi 5 (b)', 'n5-psic': 'psi 5 (c)',
    'n6-octa': 'Octahedral state', 'n6-c4': 'C4 state (pyramid)',
    'n6-omega4': 'Omega 4', 'n6-psi6': 'psi 6', 'n6-chi6': 'chi 6',
    'n6-chi7-10': 'chi 7/10 — C5 junction', 'n6-chi1-2': 'chi 1/2',
    'n6-ac2': 'AC2', 'n6-psi0': 'psi 0', 'n6-xi': 'Xi star',
    'n6-omega': 'omega', 'n6-estar': 'E star',
    'n6-eta6': 'eta 6 — sharp C2 endpoint'
  };

  function plainName(option, original) {
    const match = /^n([3-6])-(product|ghz|w)$/.exec(option.value);
    if (match) return `${{product: 'Product state', ghz: 'GHZ state', w: 'W state'}[match[2]]} (N = ${match[1]})`;
    const dicke = /^n([3-6])-dicke-?(\d+)$/.exec(option.value);
    if (dicke) return `Dicke state (N = ${dicke[1]}, k = ${dicke[2]})`;
    return plainNames[option.value] || original;
  }

  function typeset(id, fallback) {
    return window.AtlasMath ? AtlasMath.html(AtlasMath.state(id, fallback)) : fallback;
  }

  function isOpen() { return !menu.hidden; }

  function setActive(index, scroll = true) {
    if (!options.length) return;
    active = Math.max(0, Math.min(options.length - 1, index));
    optionNodes.forEach((node, i) => node.classList.toggle('is-active', i === active));
    trigger.setAttribute('aria-activedescendant', optionNodes[active].id);
    if (scroll) {
      const node = optionNodes[active], top = node.offsetTop, bottom = top + node.offsetHeight;
      if (top < menu.scrollTop) menu.scrollTop = top;
      else if (bottom > menu.scrollTop + menu.clientHeight) menu.scrollTop = bottom - menu.clientHeight;
    }
  }

  function positionMenu() {
    const rect = trigger.getBoundingClientRect();
    if (rect.bottom < 0 || rect.top > window.innerHeight) { close(); return; }
    const below = window.innerHeight - rect.bottom - 14, above = rect.top - 14;
    const upwards = below < 220 && above > below;
    wrapper.classList.toggle('opens-up', upwards);
    menu.style.maxHeight = `${Math.max(0, Math.min(420, upwards ? above : below))}px`;
  }

  function open() {
    rebuild();
    if (!options.length || select.disabled) return;
    menu.hidden = false;
    trigger.setAttribute('aria-expanded', 'true');
    positionMenu();
    setActive(Math.max(0, options.findIndex(option => option.value === select.value)));
  }

  function close() {
    menu.hidden = true;
    trigger.setAttribute('aria-expanded', 'false');
    trigger.removeAttribute('aria-activedescendant');
    search = '';
  }

  function choose(index) {
    if (!options[index]) return;
    select.value = options[index].value;
    close();
    select.dispatchEvent(new Event('change', {bubbles: true}));
    trigger.focus({preventScroll: true});
  }

  function rebuild() {
    for (const option of select.options) {
      if (!originalLabels.has(option)) originalLabels.set(option, option.textContent);
      const text = plainName(option, originalLabels.get(option));
      if (text !== option.textContent) option.textContent = text;
    }
    const next = [...select.options].filter(option => option.value && !option.disabled);
    const signature = next.map(option => `${option.value}\u0000${originalLabels.get(option)}`).join('\u0001');
    options = next;
    if (signature === optionSignature) return;
    optionSignature = signature;
    optionNodes = [];
    menu.replaceChildren();
    let previousGroup = null;
    options.forEach((option, index) => {
      const group = option.parentElement.tagName === 'OPTGROUP' ? option.parentElement.label : null;
      if (group && group !== previousGroup) {
        const heading = document.createElement('div');
        heading.className = 'state-picker-group';
        heading.setAttribute('role', 'presentation');
        heading.textContent = group;
        menu.append(heading);
      }
      previousGroup = group;
      const node = document.createElement('div');
      node.id = `state-picker-option-${index}`;
      node.className = 'state-picker-option';
      node.setAttribute('role', 'option');
      node.setAttribute('aria-label', option.textContent);
      node.dataset.value = option.value;
      node.innerHTML = typeset(/^n5-top-/.test(option.value) ? '' : option.value, originalLabels.get(option));
      node.addEventListener('pointermove', () => setActive(index, false));
      node.addEventListener('pointerdown', event => event.preventDefault());
      node.addEventListener('click', () => choose(index));
      optionNodes.push(node);
      menu.append(node);
    });
  }

  function paint() {
    const option = select.selectedOptions[0];
    const state = currentState && currentState.id === select.value ? currentState : null;
    const fallback = state?.mathLabel || (option ? originalLabels.get(option) || option.textContent : emptyText);
    if (!select.value) display.textContent = emptyText;
    else display.innerHTML = typeset(select.value, fallback);
    optionNodes.forEach((node, i) => node.setAttribute('aria-selected', String(options[i].value === select.value)));
    trigger.disabled = select.disabled;
  }

  function onKey(event) {
    if (event.key === 'Tab') { close(); return; }
    if (event.key === 'Escape') { if (isOpen()) { event.preventDefault(); close(); } return; }
    if (event.key === 'ArrowDown' || event.key === 'ArrowUp') {
      event.preventDefault();
      if (!isOpen()) open();
      else setActive(active + (event.key === 'ArrowDown' ? 1 : -1));
      return;
    }
    if (event.key === 'Home' || event.key === 'End') {
      event.preventDefault(); if (!isOpen()) open();
      setActive(event.key === 'Home' ? 0 : options.length - 1); return;
    }
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault(); if (isOpen()) choose(active); else open(); return;
    }
    if (event.key.length === 1 && !event.ctrlKey && !event.metaKey && !event.altKey) {
      event.preventDefault();
      const now = Date.now();
      search = (now - searchTime > 700 ? '' : search) + event.key.toLowerCase(); searchTime = now;
      if (!isOpen()) open();
      const query = search.normalize('NFKD');
      const start = query.length > 1 ? 0 : 1;
      for (let step = start; step < options.length + start; step++) {
        const index = (active + step) % options.length;
        if (options[index].textContent.toLowerCase().normalize('NFKD').startsWith(query)) { setActive(index); break; }
      }
    }
  }

  function init() {
    if (select) return true;
    select = document.getElementById('state-select');
    if (!select) return false;
    wrapper = document.createElement('div'); wrapper.className = 'state-picker';
    select.parentNode.insertBefore(wrapper, select); wrapper.append(select);
    trigger = document.createElement('button'); trigger.type = 'button'; trigger.id = 'state-picker-trigger';
    trigger.className = 'state-picker-trigger'; trigger.setAttribute('role', 'combobox');
    trigger.setAttribute('aria-label', 'Named state'); trigger.setAttribute('aria-haspopup', 'listbox');
    trigger.setAttribute('aria-expanded', 'false'); trigger.setAttribute('aria-controls', 'state-picker-options');
    display = document.createElement('span'); display.className = 'state-picker-display'; trigger.append(display);
    menu = document.createElement('div'); menu.id = 'state-picker-options'; menu.className = 'state-picker-options';
    menu.setAttribute('role', 'listbox'); menu.setAttribute('aria-label', 'Named states'); menu.hidden = true;
    wrapper.append(trigger, menu);
    select.hidden = true; select.tabIndex = -1; select.setAttribute('aria-hidden', 'true');
    document.querySelector('label[for="state-select"]')?.setAttribute('for', trigger.id);
    trigger.addEventListener('click', () => isOpen() ? close() : open());
    trigger.addEventListener('keydown', onKey);
    trigger.addEventListener('blur', close);
    select.addEventListener('change', () => { currentState = null; rebuild(); paint(); });
    document.addEventListener('pointerdown', event => { if (!wrapper.contains(event.target)) close(); });
    window.addEventListener('resize', () => { if (isOpen()) positionMenu(); });
    window.addEventListener('scroll', () => { if (isOpen()) positionMenu(); }, {passive: true});
    new MutationObserver(() => { rebuild(); paint(); }).observe(select, {childList: true, subtree: true});
    rebuild(); paint();
    return true;
  }

  function update(state) {
    if (!init()) return;
    currentState = state || null; rebuild(); paint();
    if (isOpen()) setActive(Math.max(0, options.findIndex(option => option.value === select.value)));
  }

  function clear(text = 'Explore a point on the plot') {
    if (!init()) return;
    currentState = null; emptyText = text; close(); paint();
  }

  window.StatePicker = {update, clear};
  init();
})();
