/* Commec Screening Report — self-contained renderer.
 *
 * This file is where most of the report is actually constructed. The Python side
 * (json_html_output.py) only flattens the ScreenResult into window.COMMEC_REPORT and
 * inlines these assets into one HTML file; everything the reader sees and interacts
 * with is built here in the browser: grouping hits per organism, computing each hit's
 * disposition, the status tabs, the sequence table, the expandable detail view
 * (swimlane lanes, hit rail, best-target card, control-list drill-down), and all
 * click handling.
 *
 * Ported from the Claude Design "DCLogic" prototype. The interactive-framework
 * pieces (x-dc / sc-for / sc-if / {{ }} / props / scope toggle) are replaced by a
 * plain string-builder + a single delegated click handler. Disposition and tab
 * status are derived from commec's own authoritative per-hit statuses (no scope
 * re-classification), so the report always agrees with commec's flag counts.
 *
 * Input: window.COMMEC_REPORT = { meta, lists, sequences: [ { ..., hits:[...] } ] }
 */
(function () {
  "use strict";

  // The logo markup is injected at HTML-assembly time from report_assets/commec-logo.svg.
  var COMMEC_LOGO = "__COMMEC_LOGO__";

  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  }

  function Report(root, data) {
    this.root = root;
    this.data = data;
    this.state = { tab: null, open: {}, selHit: {}, selList: {} };
    this._cbs = [];
  }

  Report.prototype.setState = function (patch) {
    var k;
    for (k in patch) if (patch.hasOwnProperty(k)) this.state[k] = patch[k];
    this.rerender();
  };

  // Register a click callback for this render pass; returns its index.
  Report.prototype.cb = function (fn) {
    var i = this._cbs.length;
    this._cbs.push(fn);
    return i;
  };

  /* ------------------------------------------------------------------ *
   *  Helpers (ported verbatim from the design)                          *
   * ------------------------------------------------------------------ */
  Report.prototype.num = function (n) { return (n == null || n === "") ? "" : Number(n).toLocaleString("en-US"); };
  Report.prototype.isSkip = function (s) { return (s.screenStatus || "").indexOf("Skip") === 0; };
  Report.prototype.ev = function (x) { var v = parseFloat(x); return isNaN(v) ? 1 : v; };

  // Disposition of a single hit, from commec's own status + annotations (no scope).
  //   flagbio / flagbest / warning — commec flagged (or warned on) this hit.
  //   cleared  — commec found a regulated match but cleared it as common / low-concern
  //              (rawStatus "Flag (Cleared)" / "Warning (Cleared)").
  //   exempt   — passed: matched only control lists outside the screened regions.
  //   common   — housekeeping / rRNA / common vector (low-concern search).
  Report.prototype.dispo = function (h) {
    if (h.step.indexOf("Biorisk") === 0) return (h.regulated === "Virulence Factor") ? "warning" : "flagbio";
    if (h.step.indexOf("Low Concern") === 0) return "common";
    var hasLists = h.lists && h.lists.length;
    if (h.rawStatus === "Flag") return "flagbest";
    if (h.rawStatus === "Warning") return "warning";
    if (h.rawStatus === "Flag (Cleared)" || h.rawStatus === "Warning (Cleared)") return "cleared";
    // Pass:
    if (hasLists) return "exempt";
    return "common";
  };
  Report.prototype.prio = function (d) { return ["flagbio", "warning", "flagbest", "exempt", "cleared", "common"].indexOf(d); };
  Report.prototype.dotFor = function (d) { return { flagbio: "#C23A14", warning: "#E0A020", flagbest: "#F05023", exempt: "#419BB9", cleared: "#7FA8B8", common: "#B9C0CC" }[d] || "#B9C0CC"; };
  Report.prototype.dispoLabel = function (d) { return { flagbio: "Biorisk", warning: "Virulence factor", flagbest: "Flagged", exempt: "Regulated elsewhere", cleared: "Cleared", common: "Common" }[d] || "Recorded"; };
  Report.prototype.dispoBadgeBg = function (d) {
    if (d === "flagbio" || d === "flagbest") return "#F05023";
    if (d === "warning") return "#D98A00";
    if (d === "exempt") return "#419BB9";
    if (d === "cleared") return "#7FA8B8";
    return "#8a8f9e";
  };
  Report.prototype.stepShort = function (step) { if (step.indexOf("Nucleotide") === 0) return "Nucleotide search"; if (step.indexOf("Protein") === 0) return "Protein search"; if (step.indexOf("Biorisk") === 0) return "Biorisk HMM"; return step; };
  Report.prototype.shortName = function (n) { return (n || "").replace(/\s*\(.*?\)\s*$/, "").trim(); };

  Report.prototype.geom = function (h, len) {
    var left = ((h.qs - 1) / len) * 100, w = ((h.qe - h.qs + 1) / len) * 100;
    if (isNaN(left) || left < 0) left = 0;
    if (isNaN(w) || w < 0.6) w = 0.6;
    if (left + w > 100) left = Math.max(0, 100 - w);
    return { left: left, w: w };
  };
  Report.prototype.coord = function (h) { return this.num(h.qs) + "–" + this.num(h.qe) + " bp"; };

  // Group hits by organism, choose a representative, sort by severity + e-value.
  Report.prototype.model = function (s) {
    var self = this;
    var len = s.length || 1;
    var hits = (s.hits || []).map(function (h) { return { h: h, d: self.dispo(h) }; });
    var byKey = {}, order = [];
    hits.forEach(function (x, i) {
      var key = x.h.species || x.h.name;
      if (!byKey[key]) { byKey[key] = { key: key, name: self.shortName(x.h.name), hitIdxs: [], dispo: x.d }; order.push(key); }
      byKey[key].hitIdxs.push(i);
      if (self.prio(x.d) < self.prio(byKey[key].dispo)) byKey[key].dispo = x.d;
    });
    var groups = order.map(function (k) { return byKey[k]; });
    groups.forEach(function (g) {
      var cands = g.hitIdxs.map(function (i) { return hits[i]; }).filter(function (x) { return x.d === g.dispo; });
      cands.sort(function (a, b) { return self.ev(a.h.eValue) - self.ev(b.h.eValue); });
      g.rep = cands[0].h;
      g.seg = g.hitIdxs.length;
    });
    groups.sort(function (a, b) { return (self.prio(a.dispo) - self.prio(b.dispo)) || (self.ev(a.rep.eValue) - self.ev(b.rep.eValue)); });
    var hasFlag = groups.some(function (g) { return g.dispo === "flagbio" || g.dispo === "flagbest"; });
    var hasWarn = groups.some(function (g) { return g.dispo === "warning"; });
    return { len: len, hits: hits, groups: groups, hasFlag: hasFlag, hasWarn: hasWarn };
  };

  // Sequence tab bucket — straight from commec's authoritative screen status.
  Report.prototype.status = function (s) {
    var ss = s.screenStatus || "";
    if (ss.indexOf("Skip") === 0) return "Skipped";
    if (ss === "Error" || ss === "Incomplete") return "Skipped";
    if (ss === "Flag") return "Flag";
    if (ss.indexOf("Warning") === 0 && ss.indexOf("Cleared") < 0) return "Warning";
    return "Clear";
  };

  Report.prototype.buildDetail = function (s, m) {
    var self = this;
    var sid = s.id;
    var hits = m.hits, len = m.len;
    if (!hits.length) {
      return {
        axisMid: this.num(Math.round(len / 2)), axisEnd: this.num(len),
        bioBars: [], warnBars: [], bestBars: [], annBars: [], bioHas: false, warnHas: false, bioneither: true, bestEmpty: true, annEmpty: true,
        railSections: [], railLabel: "All hits · 0", hasHits: false, noHits: true,
        selBadge: "Clear", selBadgeBg: "#419BB9", isBiorisk: false, descMode: false, showBest: true,
        selDesc: "", selDomain: "—", selSearch: "", selName: "No matches",
        selNameStyle: "", selCategory: "—", selTaxid: "—", selTarget: "—", selTargetDesc: "",
        selPct: "—", selEval: "—", selCoords: "—",
        controlLabel: "No regions of concern", tags: [],
        drill: { show: false, tag: "", name: "", region: "", code: "", authority: "", entry: "" }
      };
    }
    var sel = this.state.selHit[sid];
    if (sel == null || sel >= hits.length) sel = hits.findIndex(function (x) { return x.d === "flagbio" || x.d === "flagbest" || x.d === "warning"; });
    if (sel < 0) sel = 0;
    var selList = (this.state.selList[sid] == null) ? -1 : this.state.selList[sid];

    var selectHit = function (i) {
      var sh = {}, sl = {}; var k;
      for (k in self.state.selHit) sh[k] = self.state.selHit[k];
      for (k in self.state.selList) sl[k] = self.state.selList[k];
      sh[sid] = i; sl[sid] = -1;
      self.setState({ selHit: sh, selList: sl });
    };

    // lanes — one bar per hit, selection keyed to hit index
    var bio = [], warn = [], best = [], ann = [];
    hits.forEach(function (x, i) {
      var on = i === sel;
      var g = self.geom(x.h, len);
      var ring = on ? "box-shadow:0 0 0 2px #23285A;" : "";
      var mk = function (top, h, bg, zBoost) {
        var zi = on ? 9999 : Math.max(2, Math.round(300 - g.w * 2)) + (zBoost || 0);
        return { style: "position:absolute; top:" + top + "px; height:" + h + "px; left:" + g.left + "%; width:" + g.w + "%; min-width:6px; background:" + bg + "; opacity:.82; border-radius:2px; z-index:" + zi + "; cursor:pointer; " + ring, title: self.shortName(x.h.name) + " · " + (x.h.pctId || "") + " · " + self.coord(x.h), onSelect: function () { selectHit(i); } };
      };
      if (x.d === "flagbio") bio.push(mk(6, 22, "#C23A14", 500));
      else if (x.d === "warning") warn.push(mk(6, 22, "#E0A020", 500));
      else if (x.d === "flagbest") { var b = mk(9, 22, "#F05023"); b.label = self.shortName(x.h.name); b.labelStyle = "position:absolute; top:12px; left:calc(" + g.left + "% + " + (g.w > 6 ? 8 : 9) + "px); font-family:Manrope,Arial,sans-serif; text-transform:uppercase; letter-spacing:.04em; font-size:8px; font-weight:700; color:#C23A14; white-space:nowrap; pointer-events:none;"; best.push(b); }
      else if (x.d === "exempt") ann.push(mk(9, 13, "rgba(65,155,185,.6)"));
      else if (x.d === "cleared") ann.push(mk(9, 13, "rgba(127,168,184,.45)"));
      else ann.push(mk(9, 13, "rgba(185,192,204,.85)"));
    });

    // rail — one row per hit, split into Flag / Warning / Clear sections
    var mkRow = function (x, i) {
      var on = i === sel, h = x.h;
      var ital = h.category === "Viruses" || h.category === "Bacteria";
      return {
        name: self.shortName(h.name) || h.name, coords: self.coord(h) + " · " + self.stepShort(h.step).replace(/ search$/i, ""), pct: h.pctId || "", dot: self.dotFor(x.d),
        weight: on ? 700 : 400, nameStyle: ital ? "font-style:italic;" : "",
        rowStyle: "display:flex; align-items:flex-start; gap:9px; padding:9px 11px; border-radius:3px; cursor:pointer; border:1px solid " + (on ? "#C7CBD8" : "transparent") + "; background:" + (on ? "#fff" : "transparent") + ";",
        onSelect: function () { selectHit(i); }
      };
    };
    var railCats = [
      { label: "Flag", color: "#F05023", match: function (dd) { return dd === "flagbio" || dd === "flagbest"; } },
      { label: "Warning", color: "#D98A00", match: function (dd) { return dd === "warning"; } },
      { label: "Clear", color: "#419BB9", match: function (dd) { return !(dd === "flagbio" || dd === "flagbest" || dd === "warning"); } }
    ];
    var railSections = railCats.map(function (c) {
      var items = hits.map(function (x, i) { return { x: x, i: i }; }).filter(function (o) { return c.match(o.x.d); }).map(function (o) { return mkRow(o.x, o.i); });
      return { label: c.label, color: c.color, count: items.length, items: items };
    }).filter(function (sec) { return sec.items.length; });

    var x = hits[sel] || hits[0];
    var h = x.h, d = x.d;
    var isBiorisk = (h.step || "").indexOf("Biorisk") === 0;
    var isLowConcern = (h.step || "").indexOf("Low Concern") === 0;
    var descMode = isBiorisk || isLowConcern;
    var allLists = h.lists || [];
    var isFlag = d === "flagbio" || d === "flagbest" || d === "warning";
    var tagLists = allLists;
    var tags = (tagLists.length ? tagLists : [{ acronym: "—", short: "—", name: "Shared protein · not diagnostic" }]).map(function (l, i) {
      var active = (i === selList) && tagLists.length > 0;
      return {
        tag: l.acronym || l.short || "—", name: self.shortName(l.name), tagColor: active ? "#23285A" : "#8a8f9e",
        rowStyle: "display:inline-flex; align-items:center; gap:6px; padding:5px 10px; border-radius:3px; cursor:pointer; border:1px solid " + (active ? "#23285A" : "#E4E7EE") + "; background:" + (active ? "#F4F6FB" : "#fff") + ";",
        onSelect: function () {
          var sl = {}; var k; for (k in self.state.selList) sl[k] = self.state.selList[k];
          sl[sid] = (selList === i ? -1 : i);
          self.setState({ selList: sl });
        }
      };
    });
    var dl = (tagLists.length && selList >= 0 && selList < tagLists.length) ? tagLists[selList] : null;
    var drill = dl
      ? { show: true, tag: dl.acronym || dl.short || "—", name: dl.acronym || dl.short || "—", region: dl.region || "—", code: dl.acronym || "—", authority: dl.name, entry: dl.source || "—" }
      : { show: false, tag: "", name: "", region: "", code: "", authority: "", entry: "" };

    var controlLabel, badge = this.dispoLabel(d);
    var listCount = function () { return allLists.length + " list" + (allLists.length > 1 ? "s" : ""); };
    if (d === "flagbest" || d === "flagbio") { controlLabel = allLists.length ? ("Controlled by " + listCount() + " — click one") : "In scope"; }
    else if (d === "warning") { controlLabel = "Virulence factor"; }
    else if (d === "exempt") { controlLabel = "Controlled by " + listCount() + " — outside the screened regions"; }
    else if (d === "cleared") { controlLabel = allLists.length ? ("Cleared as common / low-concern · matched " + listCount()) : "Cleared as common / low-concern"; }
    else { controlLabel = "Common / low-concern"; }

    return {
      axisMid: this.num(Math.round(len / 2)), axisEnd: this.num(len),
      bioBars: bio, warnBars: warn, bestBars: best, annBars: ann,
      bioHas: bio.length > 0, warnHas: warn.length > 0, bioneither: bio.length === 0 && warn.length === 0, bestEmpty: best.length === 0, annEmpty: ann.length === 0,
      railSections: railSections, railLabel: "All hits · " + hits.length, hasHits: true, noHits: false,
      selBadge: badge, selBadgeBg: this.dispoBadgeBg(d),
      isBiorisk: isBiorisk, descMode: descMode, showBest: !descMode,
      selDesc: h.desc || h.description || "(no description provided)", selDomain: h.domain || "—",
      selSearch: this.stepShort(h.step), selName: this.shortName(h.name) || h.name,
      selNameStyle: (h.category === "Viruses" || h.category === "Bacteria") ? "font-style:italic;" : "",
      selCategory: h.category || "—", selTaxid: h.taxid || "—",
      selTarget: h.target || "—", selTargetDesc: h.targetDesc || "", selPct: h.pctId || "—",
      selEval: h.eValue || "—", selCoords: this.coord(h),
      controlLabel: controlLabel, tags: tags, drill: drill
    };
  };

  Report.prototype.renderVals = function () {
    var self = this;
    var ds = this.data;
    var m = ds.meta || {};

    var seqs = ds.sequences || [];
    var tagged = seqs.map(function (s) { return { s: s, status: self.status(s) }; });
    var counts = { Flag: 0, Warning: 0, Clear: 0, Skipped: 0 };
    tagged.forEach(function (t) { counts[t.status] = (counts[t.status] || 0) + 1; });

    var order = ["Flag", "Warning", "Clear", "Skipped"];
    var tab = this.state.tab;
    if (!tab) tab = order.find(function (k) { return counts[k] > 0; }) || "Flag";

    var colorFor = { Flag: "#F05023", Warning: "#D98A00", Clear: "#419BB9", Skipped: "#8a8f9e" };
    var tabs = order.map(function (k) {
      var active = k === tab, c = counts[k] || 0;
      return {
        label: k, count: c,
        onClick: function () { self.setState({ tab: k, open: {}, selHit: {}, selList: {} }); },
        style: "display:flex; align-items:baseline; gap:9px; padding:14px 20px 14px; cursor:pointer; border-bottom:2px solid " + (active ? "#23285A" : "transparent") + "; background:" + (active ? "#F8F9FC" : "transparent") + ";",
        numStyle: "font-family:'Crimson Pro',Georgia,serif; font-size:30px; font-weight:700; line-height:1; color:" + colorFor[k] + ";" + (c === 0 && !active ? "opacity:.32;" : ""),
        labelStyle: "font-size:11px; font-weight:700; letter-spacing:.1em; color:#23285A;" + (active ? "" : "opacity:.6;")
      };
    });

    var filtered = tagged.filter(function (t) { return t.status === tab; });
    var rows = filtered.map(function (t) {
      var s = t.s, sid = s.id;
      var skip = t.status === "Skipped";
      var open = !!self.state.open[sid] && !skip;
      var lengthLabel = self.num(s.length) + " bp";

      var finding = "", mapBars = [];
      if (!skip) {
        var mm = self.model(s);
        var isBio = function (x) { return (x.h.step || "").indexOf("Biorisk") === 0; };
        var buildTrack = function (sevOf, colorOf, top, hgt) {
          var N = 120, bins = new Array(N).fill(0);
          mm.hits.forEach(function (x) {
            var sev = sevOf(x); if (!sev) return;
            var b0 = Math.floor(Math.max(0, (x.h.qs - 1) / mm.len) * N);
            var b1 = Math.ceil(Math.min(1, x.h.qe / mm.len) * N);
            if (b1 <= b0) b1 = b0 + 1;
            for (var b = b0; b < b1 && b < N; b++) if (sev > bins[b]) bins[b] = sev;
          });
          for (var b = 0; b < N;) {
            var sev = bins[b];
            if (!sev) { b++; continue; }
            var e = b; while (e < N && bins[e] === sev) e++;
            var left = (b / N) * 100, w = ((e - b) / N) * 100;
            mapBars.push({ style: "position:absolute; top:" + top + "px; height:" + hgt + "px; left:" + left + "%; width:" + w + "%; background:" + colorOf(sev) + "; border-radius:1px;" });
            b = e;
          }
        };
        buildTrack(function (x) { return isBio(x) ? 0 : (x.d === "flagbest" ? 2 : 1); }, function (sev) { return sev === 2 ? "#F05023" : "rgba(65,155,185,.5)"; }, 2, 6);
        buildTrack(function (x) { return !isBio(x) ? 0 : (x.d === "flagbio" ? 2 : 1); }, function (sev) { return sev === 2 ? "#C23A14" : "#E0A020"; }, 10, 3);
        if (t.status === "Flag") { var fg = mm.groups.filter(function (g) { return g.dispo === "flagbest" || g.dispo === "flagbio"; }).map(function (g) { return g.dispo === "flagbio" ? "Biorisk flag" : g.name; }); var uniq = []; fg.forEach(function (n) { if (uniq.indexOf(n) < 0) uniq.push(n); }); finding = uniq.slice(0, 3).join(", ") + (uniq.length > 3 ? " +" + (uniq.length - 3) : ""); }
        else if (t.status === "Warning") finding = mm.hasWarn ? "Virulence factor · review" : "No matches in scope · risk unknown";
        else if (!mm.hits.length) { finding = s.rationale || "No hits"; }
        else {
          var count = function (dd) { return mm.groups.filter(function (g) { return g.dispo === dd; }).length; };
          var cl = count("cleared"), ex = count("exempt");
          var bits = [];
          if (cl) bits.push(cl + " cleared");
          if (ex) bits.push(ex + " regulated elsewhere");
          finding = bits.length ? bits.join(" · ") : (s.rationale || "Nothing regulated in scope");
        }
      } else { finding = s.rationale || "Not screened"; }

      var detail = open ? self.buildDetail(s, self.model(s)) : null;

      return {
        name: s.name, lengthLabel: lengthLabel, finding: finding, mapBars: mapBars, open: open,
        wrapBg: open ? "#FBFCFE" : "#fff",
        rowStyle: "display:grid; grid-template-columns:18px 3.1fr 2.4fr 3.0fr; align-items:center; gap:18px; padding:13px 40px; " + (skip ? "cursor:default;" : "cursor:pointer;") + (open ? "background:#FBFCFE;" : ""),
        chevStyle: "font-size:12px; color:" + (open ? "#F05023" : "#9a9fb0") + "; transition:transform .18s cubic-bezier(.2,.7,.2,1); transform:rotate(" + (open ? 90 : 0) + "deg);" + (skip ? "visibility:hidden;" : ""),
        onToggle: skip ? function () {} : (function () {
          var o = {}; var k; for (k in self.state.open) o[k] = self.state.open[k];
          o[sid] = !self.state.open[sid];
          self.setState({ open: o });
        }),
        detail: detail
      };
    });

    return {
      mFile: m.file || "Screening report",
      mLine1: (m.nQueries != null ? m.nQueries : seqs.length) + " sequences · " + this.num(m.totalLen) + " bp · screened " + (m.date || "") + " · run " + (m.time || ""),
      mLine2: "commec " + (m.version || "") + " · output " + (m.output || "") + (m.tools ? " · " + m.tools : ""),
      tabs: tabs, rows: rows,
      hasRows: rows.length > 0, isEmpty: rows.length === 0,
      emptyText: "No sequences with " + tab.toLowerCase() + " status in this run.",
      listCount: (ds.lists || []).length
    };
  };

  /* ------------------------------------------------------------------ *
   *  Template (string builder replacing x-dc / sc-for / sc-if)          *
   * ------------------------------------------------------------------ */
  Report.prototype.template = function (v) {
    var self = this;

    var masthead =
      '<div style="background:#23285A; padding:26px 40px 24px;">' +
        '<div style="display:flex; align-items:center; gap:18px;">' + COMMEC_LOGO +
          '<div style="flex:1;"></div>' +
          '<span class="cap" style="font-size:9.5px; font-weight:700; letter-spacing:.12em; color:rgba(255,255,255,.5);">Synthesis screening report</span>' +
        '</div>' +
        '<div style="display:flex; align-items:flex-end; justify-content:space-between; margin-top:22px; gap:30px;">' +
          '<div>' +
            '<div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:32px; font-weight:700; color:#fff; line-height:1.04; letter-spacing:-.01em;">' + esc(v.mFile) + '</div>' +
            '<div class="mono" style="font-size:11.5px; color:rgba(255,255,255,.62); margin-top:11px; line-height:1.75;">' + esc(v.mLine1) + '<br>' + esc(v.mLine2) + '</div>' +
          '</div>' +
        '</div>' +
      '</div>';

    var tabsRow =
      '<div style="display:flex; align-items:stretch; gap:0; padding:0 40px; border-bottom:1px solid #E4E7EE; background:#fff;">' +
        v.tabs.map(function (tb) {
          return '<div class="tab-hover" data-cb="' + self.cb(tb.onClick) + '" style="' + tb.style + '">' +
            '<span style="' + tb.numStyle + '">' + esc(tb.count) + '</span>' +
            '<span class="cap" style="' + tb.labelStyle + '">' + esc(tb.label) + '</span>' +
          '</div>';
        }).join("") +
        '<div style="flex:1;"></div>' +
        '<div style="display:flex; flex-direction:column; align-items:flex-end; justify-content:center; gap:3px; padding:12px 0;">' +
          '<span class="cap" style="font-size:9px; font-weight:700; letter-spacing:.1em; color:#9a9fb0;">Control lists</span>' +
          '<span class="sans" style="font-size:12px; font-weight:700; color:#23285A;">' + esc(v.listCount) + ' screened</span>' +
        '</div>' +
      '</div>';

    var tableHeader = v.hasRows ?
      '<div class="cap" style="display:grid; grid-template-columns:18px 3.1fr 2.4fr 3.0fr; align-items:center; gap:18px; padding:12px 40px 10px; font-size:10px; font-weight:700; letter-spacing:.1em; color:#9a9fb0; border-bottom:1px solid #E4E7EE;">' +
        '<div></div><div>Sequence</div><div>Map</div><div>Finding</div></div>' : "";

    var rows = v.rows.map(function (r) { return self.rowTemplate(r); }).join("");

    var empty = v.isEmpty ?
      '<div style="padding:46px 40px; text-align:center;"><div class="sans" style="font-size:13px; color:#9a9fb0;">' + esc(v.emptyText) + '</div></div>' : "";

    return '<div class="cm report-card"><div>' + masthead + tabsRow + tableHeader + rows + empty + this.howToRead() + '</div></div>';
  };

  Report.prototype.rowTemplate = function (r) {
    var self = this;
    var mapBars = r.mapBars.map(function (b) { return '<div style="' + b.style + '"></div>'; }).join("");
    var head =
      '<div class="row-hover" data-cb="' + self.cb(r.onToggle) + '" style="' + r.rowStyle + '">' +
        '<span style="' + r.chevStyle + '">▸</span>' +
        '<div style="min-width:0;">' +
          '<div style="font-family:\'Crimson Pro\',Georgia,serif; font-weight:700; font-size:14.5px; white-space:nowrap; overflow:hidden; text-overflow:ellipsis;">' + esc(r.name) + '</div>' +
          '<div class="mono" style="font-size:10.5px; color:#9a9fb0; margin-top:2px;">' + esc(r.lengthLabel) + '</div>' +
        '</div>' +
        '<div style="position:relative; height:14px; background:#F4F5FA; border:1px solid #E7EAF1; border-radius:2px;">' + mapBars + '</div>' +
        '<div class="sans" style="font-size:12px; color:#5b6070; line-height:1.4; overflow:hidden; text-overflow:ellipsis;">' + esc(r.finding) + '</div>' +
      '</div>';

    var expanded = (r.open && r.detail) ? this.detailTemplate(r.detail) : "";
    return '<div style="border-bottom:1px solid #EEF0F5; background:' + r.wrapBg + ';">' + head + expanded + '</div>';
  };

  Report.prototype.detailTemplate = function (d) {
    var self = this;

    var bioLane = d.bioHas
      ? '<div style="position:relative; min-height:30px; display:flex; align-items:center;">' + d.bioBars.map(function (b) { return '<div title="' + esc(b.title) + '" data-cb="' + self.cb(b.onSelect) + '" style="' + b.style + '"></div>'; }).join("") + '</div>'
      : "";
    var warnLane = d.warnHas
      ? '<div style="position:relative; min-height:30px; display:flex; align-items:center;">' + d.warnBars.map(function (b) { return '<div title="' + esc(b.title) + '" data-cb="' + self.cb(b.onSelect) + '" style="' + b.style + '"></div>'; }).join("") + '</div>'
      : "";
    var neitherLane = d.bioneither
      ? '<div style="position:relative; min-height:30px; display:flex; align-items:center;"><div style="position:absolute; left:0; right:0; top:50%; height:1px; background:repeating-linear-gradient(90deg,#E0E2EA,#E0E2EA 4px,transparent 4px,transparent 8px);"></div><span class="sans" style="position:relative; font-size:10px; color:#aab; background:#FBFCFE; padding-right:8px;">No biorisk models matched</span></div>'
      : "";

    var bestBars = d.bestBars.map(function (b) {
      return '<div title="' + esc(b.title) + '" data-cb="' + self.cb(b.onSelect) + '" style="' + b.style + '"></div>';
    }).join("");
    var bestEmpty = d.bestEmpty ? '<div style="position:absolute; left:0; right:0; top:50%; height:1px; background:repeating-linear-gradient(90deg,#E0E2EA,#E0E2EA 4px,transparent 4px,transparent 8px);"></div><span class="sans" style="position:absolute; top:50%; transform:translateY(-50%); font-size:10px; color:#aab; background:#FBFCFE; padding-right:8px;">Nothing flagged in scope</span>' : "";

    var annBars = d.annBars.map(function (b) { return '<div title="' + esc(b.title) + '" data-cb="' + self.cb(b.onSelect) + '" style="' + b.style + '"></div>'; }).join("");
    var annEmpty = d.annEmpty ? '<span class="sans" style="position:relative; font-size:10px; color:#aab;">None</span>' : "";

    var swimlane =
      '<div style="padding:16px 40px 14px;">' +
        '<div style="position:relative; height:13px; margin-left:132px; border-bottom:1px solid #D6DAE0;">' +
          '<span class="mono" style="position:absolute; left:0; bottom:2px; font-size:9px; color:#9a9fb0;">1</span>' +
          '<span class="mono" style="position:absolute; left:50%; bottom:2px; font-size:9px; color:#9a9fb0; transform:translateX(-50%);">' + esc(d.axisMid) + '</span>' +
          '<span class="mono" style="position:absolute; right:0; bottom:2px; font-size:9px; color:#9a9fb0;">' + esc(d.axisEnd) + '</span>' +
        '</div>' +
        '<div style="display:flex; align-items:stretch; min-height:36px; border-bottom:1px solid #F2F3F7;">' +
          '<div style="width:132px; flex-shrink:0; display:flex; flex-direction:column; justify-content:center; gap:2px;"><div style="display:flex; align-items:center; gap:7px;"><span style="width:9px;height:9px;border-radius:2px;background:#C23A14;"></span><span class="cap" style="font-size:9.5px; font-weight:700; color:#23285A;">Biorisk</span></div></div>' +
          '<div style="flex:1; display:flex; flex-direction:column; justify-content:center;">' + bioLane + warnLane + neitherLane + '</div>' +
        '</div>' +
        '<div style="display:flex; align-items:stretch; min-height:42px; border-bottom:1px solid #F2F3F7;">' +
          '<div style="width:132px; flex-shrink:0; display:flex; flex-direction:column; justify-content:center; gap:2px;"><div style="display:flex; align-items:center; gap:7px;"><span style="width:9px;height:9px;border-radius:2px;background:#F05023;"></span><span class="cap" style="font-size:9.5px; font-weight:700; color:#23285A;">Best match</span></div><span class="sans" style="font-size:8.5px; color:#9a9fb0; padding-left:16px;">click a hit</span></div>' +
          '<div style="position:relative; flex:1; min-height:42px;">' + bestEmpty + bestBars + '</div>' +
        '</div>' +
        '<div style="display:flex; align-items:stretch; min-height:40px;">' +
          '<div style="width:132px; flex-shrink:0; display:flex; flex-direction:column; justify-content:center; gap:2px;"><div style="display:flex; align-items:center; gap:7px;"><span style="width:9px;height:9px;border-radius:2px;background:#419BB9;"></span><span class="cap" style="font-size:9.5px; font-weight:700; color:#23285A;">Annotations</span></div><span class="sans" style="font-size:8.5px; color:#9a9fb0; padding-left:16px;">context</span></div>' +
          '<div style="position:relative; flex:1; min-height:40px; display:flex; align-items:center;">' + annEmpty + annBars + '</div>' +
        '</div>' +
      '</div>';

    var body = "";
    if (d.hasHits) body = this.detailBody(d);
    else if (d.noHits) body = '<div style="border-top:1px solid #E4E7EE; padding:20px 40px; background:#F8F9FC;"><div class="sans" style="font-size:12px; color:#8a8f9e; font-style:italic;">No hits — nothing matched any screening database for this sequence.</div></div>';

    return '<div style="background:#FBFCFE; border-top:1px solid #E7EAF1;">' + swimlane + body + '</div>';
  };

  Report.prototype.detailBody = function (d) {
    var self = this;

    var rail =
      '<div style="border-right:1px solid #E4E7EE; background:#F8F9FC; padding:12px 10px;">' +
        '<div class="cap" style="font-size:9px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; padding:2px 4px 9px;">' + esc(d.railLabel) + '</div>' +
        d.railSections.map(function (sec) {
          return '<div style="display:flex; align-items:center; gap:6px; padding:11px 4px 5px;">' +
              '<span style="width:7px;height:7px;border-radius:2px;background:' + sec.color + ';flex-shrink:0;"></span>' +
              '<span class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:' + sec.color + ';">' + esc(sec.label) + '</span>' +
              '<span class="mono" style="font-size:9px; color:#b3b8c4;">' + esc(sec.count) + '</span>' +
            '</div>' +
            '<div style="display:flex; flex-direction:column; gap:2px;">' +
              sec.items.map(function (rl) {
                return '<div data-cb="' + self.cb(rl.onSelect) + '" style="' + rl.rowStyle + '">' +
                  '<span style="width:8px;height:8px;border-radius:9999px;background:' + rl.dot + ';flex-shrink:0;margin-top:4px;"></span>' +
                  '<div style="flex:1; min-width:0;"><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; font-weight:' + rl.weight + '; ' + rl.nameStyle + ' line-height:1.25;">' + esc(rl.name) + '</div><div class="mono" style="font-size:9.5px; color:#9a9fb0; margin-top:1px;">' + esc(rl.coords) + '</div></div>' +
                  '<span class="mono" style="font-size:10.5px; color:#6b7080;">' + esc(rl.pct) + '</span>' +
                '</div>';
              }).join("") +
            '</div>';
        }).join("") +
      '</div>';

    var meta = "";
    if (d.showBest) {
      meta = '<div class="sans" style="font-size:11.5px; color:#6b7080; margin-top:5px; display:flex; gap:14px; flex-wrap:wrap;">' +
        '<span>Category <strong style="color:#23285A;">' + esc(d.selCategory) + '</strong></span>' +
        '<span>taxid <span class="mono" style="color:#23285A;">' + esc(d.selTaxid) + '</span></span>' +
      '</div>';
    } else if (d.isBiorisk) {
      meta = '<div class="sans" style="font-size:11.5px; color:#6b7080; margin-top:5px; display:flex; gap:14px; flex-wrap:wrap;">' +
        '<span>Domain <strong style="color:#23285A;">' + esc(d.selDomain) + '</strong></span></div>';
    }

    var hero = "";
    if (d.showBest) {
      hero = '<div style="margin:16px 0; padding:20px; background:#fff; border:1px solid #E4E7EE; border-radius:4px; text-align:center;">' +
        '<div class="cap" style="font-size:9px; font-weight:700; letter-spacing:.12em; color:#8a8f9e;">Best target</div>' +
        '<div class="mono" style="font-size:23px; font-weight:600; color:#419BB9; margin-top:8px; letter-spacing:.01em;">' + esc(d.selTarget) + '</div>' +
        '<div class="sans" style="font-size:12.5px; color:#3c4156; margin-top:6px; max-width:520px; margin-left:auto; margin-right:auto; line-height:1.45;">' + esc(d.selTargetDesc) + '</div>' +
        '<div style="display:flex; justify-content:center; gap:34px; margin-top:16px; padding-top:15px; border-top:1px solid #E7EAF1;">' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; margin-bottom:5px;">Identity</div><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:20px; font-weight:700; color:#23285A;">' + esc(d.selPct) + '</div></div>' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; margin-bottom:5px;">E-value</div><div class="mono" style="font-size:16px; font-weight:600; color:#23285A; margin-top:2px;">' + esc(d.selEval) + '</div></div>' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; margin-bottom:5px;">Region</div><div class="mono" style="font-size:13px; color:#23285A; margin-top:5px;">' + esc(d.selCoords) + '</div></div>' +
        '</div>' +
      '</div>';
    } else if (d.descMode) {
      hero = '<div style="margin:16px 0; padding:20px; background:#fff; border:1px solid #E4E7EE; border-radius:4px;">' +
        '<div class="cap" style="font-size:9px; font-weight:700; letter-spacing:.12em; color:#8a8f9e;">Description</div>' +
        '<div class="sans" style="font-size:13.5px; color:#3c4156; margin-top:8px; line-height:1.5;">' + esc(d.selDesc) + '</div>' +
        '<div style="display:flex; gap:44px; margin-top:16px; padding-top:15px; border-top:1px solid #E7EAF1;">' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; margin-bottom:5px;">E-value</div><div class="mono" style="font-size:16px; font-weight:600; color:#23285A; margin-top:2px;">' + esc(d.selEval) + '</div></div>' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; margin-bottom:5px;">Region</div><div class="mono" style="font-size:13px; color:#23285A; margin-top:5px;">' + esc(d.selCoords) + '</div></div>' +
        '</div>' +
      '</div>';
    }

    var tagsBlock = "";
    if (d.showBest) {
      tagsBlock = '<div class="cap" style="font-size:9px; font-weight:700; letter-spacing:.1em; color:#8a8f9e; margin-bottom:9px;">' + esc(d.controlLabel) + '</div>' +
        '<div style="display:flex; flex-wrap:wrap; gap:6px;">' +
          d.tags.map(function (t) {
            return '<div data-cb="' + self.cb(t.onSelect) + '" style="' + t.rowStyle + '">' +
              '<span class="cap" style="font-size:9px; font-weight:700; letter-spacing:.05em; color:' + t.tagColor + ';">' + esc(t.tag) + '</span>' +
              '<span class="sans" style="font-size:11px; color:#6b7080;">' + esc(t.name) + '</span>' +
            '</div>';
          }).join("") +
        '</div>';
    }

    var drill = "";
    if (d.drill.show) {
      drill = '<div style="margin-top:13px; border:1px solid #D7DCE8; border-radius:4px; overflow:hidden;">' +
        '<div style="display:flex; align-items:center; gap:9px; padding:11px 14px; background:#23285A;">' +
          '<span class="cap" style="font-size:9.5px; font-weight:700; letter-spacing:.05em; color:#23285A; background:#fff; padding:3px 7px; border-radius:2px;">' + esc(d.drill.tag) + '</span>' +
          '<span class="sans" style="font-size:12.5px; font-weight:600; color:#fff;">' + esc(d.drill.name) + '</span>' +
        '</div>' +
        '<div style="padding:13px 15px; display:grid; grid-template-columns:1fr 1fr; gap:12px 20px; background:#fff;">' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#9a9fb0; margin-bottom:3px;">Jurisdiction</div><div class="sans" style="font-size:12px; color:#23285A;">' + esc(d.drill.region) + '</div></div>' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#9a9fb0; margin-bottom:3px;">Control list</div><div class="sans" style="font-size:12px; color:#23285A; line-height:1.4;">' + esc(d.drill.authority) + '</div></div>' +
          '<div><div class="cap" style="font-size:8.5px; font-weight:700; letter-spacing:.1em; color:#9a9fb0; margin-bottom:3px;">Listed entry</div><div class="sans" style="font-size:12px; color:#23285A; font-style:italic;">' + esc(d.drill.entry) + '</div></div>' +
        '</div>' +
      '</div>';
    }

    var pane =
      '<div style="padding:20px 26px 24px;">' +
        '<div style="display:flex; align-items:center; gap:9px; margin-bottom:9px;">' +
          '<span class="badge" style="background:' + d.selBadgeBg + '; color:#fff;">' + esc(d.selBadge) + '</span>' +
          '<span class="sans" style="font-size:11px; color:#6b7080;">' + esc(d.selSearch) + '</span>' +
        '</div>' +
        '<div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:21px; font-weight:700; ' + d.selNameStyle + ' line-height:1.12;">' + esc(d.selName) + '</div>' +
        meta + hero + tagsBlock + drill +
      '</div>';

    return '<div style="display:grid; grid-template-columns:250px 1fr; border-top:1px solid #E4E7EE; align-items:stretch;">' + rail + pane + '</div>';
  };

  Report.prototype.howToRead = function () {
    return '<div style="padding:24px 40px 28px; background:#F8F9FC; border-top:1px solid #E4E7EE;">' +
      '<div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:17px; font-weight:700; margin-bottom:10px;">How to read this report</div>' +
      '<p style="font-family:\'Crimson Pro\',Georgia,serif; font-size:13.5px; line-height:1.6; color:#4c5262; margin:0 0 20px; max-width:840px;">Pick a status tab, then click any sequence to see its screening results, including matching HMM models or sequences, and the control list entries behind each hit.</p>' +
      '<div style="display:grid; grid-template-columns:1fr 1fr 1fr; gap:26px;">' +
        '<div>' +
          '<div class="cap" style="font-size:10px; font-weight:700; letter-spacing:.1em; color:#23285A; padding-bottom:7px; border-bottom:1px solid #E4E7EE; margin-bottom:12px;">Status</div>' +
          '<div style="display:flex; flex-direction:column; gap:11px;">' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span class="badge" style="background:#F05023; color:#fff; margin-top:1px;">Flag</span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;">Found biorisk or match to a controlled organism. Follow-up recommended.</div></div>' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span class="badge" style="background:#D98A00; color:#fff; margin-top:1px;">Warning</span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;">Some signals of concern &mdash; consider follow-up screening.</div></div>' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span class="badge" style="background:#419BB9; color:#fff; margin-top:1px;">Clear</span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;">Nothing of concern identified in-scope; any matches were cleared.</div></div>' +
          '</div>' +
        '</div>' +
        '<div>' +
          '<div class="cap" style="font-size:10px; font-weight:700; letter-spacing:.1em; color:#23285A; padding-bottom:7px; border-bottom:1px solid #E4E7EE; margin-bottom:12px;">Screening lanes</div>' +
          '<div style="display:flex; flex-direction:column; gap:11px;">' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span style="width:11px;height:11px;border-radius:2px;background:#C23A14;margin-top:3px;flex-shrink:0;"></span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;"><strong>Biorisk</strong> — HMM screen for curated sequences of concern and virulence factors.</div></div>' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span style="width:11px;height:11px;border-radius:2px;background:#F05023;margin-top:3px;flex-shrink:0;"></span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;"><strong>Best match</strong> — closest reference is an organism controlled in the selected regions.</div></div>' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span style="width:11px;height:11px;border-radius:2px;background:#419BB9;margin-top:3px;flex-shrink:0;opacity:.7;"></span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;"><strong>Annotations</strong> — regulated-elsewhere, cleared and low-concern context.</div></div>' +
          '</div>' +
        '</div>' +
        '<div>' +
          '<div class="cap" style="font-size:10px; font-weight:700; letter-spacing:.1em; color:#23285A; padding-bottom:7px; border-bottom:1px solid #E4E7EE; margin-bottom:12px;">Annotation types</div>' +
          '<div style="display:flex; flex-direction:column; gap:11px;">' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span style="width:18px;height:11px;border-radius:2px;background:#419BB9;opacity:.55;margin-top:3px;flex-shrink:0;"></span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;"><strong>Regulated elsewhere</strong> — controlled only under lists outside the screened regions.</div></div>' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span style="width:18px;height:11px;border-radius:2px;background:rgba(127,168,184,.55);margin-top:3px;flex-shrink:0;"></span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;"><strong>Cleared</strong> — matched a controlled organism but cleared as common / low-concern.</div></div>' +
            '<div style="display:flex; gap:11px; align-items:flex-start;"><span style="width:18px;height:11px;border-radius:2px;background:#B9C0CC;margin-top:3px;flex-shrink:0;"></span><div style="font-family:\'Crimson Pro\',Georgia,serif; font-size:12.5px; line-height:1.5; color:#3c4156;"><strong>Common / low-concern</strong> — housekeeping genes, rRNA, common vectors.</div></div>' +
          '</div>' +
        '</div>' +
      '</div>' +
    '</div>';
  };

  Report.prototype.rerender = function () {
    this._cbs = [];
    var vals = this.renderVals();
    this.root.innerHTML = this.template(vals);
  };

  Report.prototype.mount = function () {
    var self = this;
    this.root.addEventListener("click", function (e) {
      var el = e.target;
      while (el && el !== self.root) {
        if (el.getAttribute && el.getAttribute("data-cb") != null) {
          var fn = self._cbs[+el.getAttribute("data-cb")];
          if (fn) fn();
          return;
        }
        el = el.parentNode;
      }
    });
    this.rerender();
  };

  function boot() {
    var root = document.getElementById("commec-report");
    if (!root) return;
    var data = window.COMMEC_REPORT || { meta: {}, lists: [], sequences: [] };
    new Report(root, data).mount();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
  else boot();
})();
