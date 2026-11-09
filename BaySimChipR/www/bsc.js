/* BaySimChip (R Shiny) - helpers for the drag-and-drop model editor.
 *
 * The model diagram is a vis.js network rendered by visNetwork; its instance
 * is available as document.getElementById('graph' + outputId).chart.
 * The toolbar buttons of the "Model structure" tab switch the network into
 * vis.js manipulation modes (add node / add edge / delete).  The add
 * callbacks (defined in R, utils_ui.R) do not add anything themselves: they
 * send the request to the server as input$<outputId>_graphChange, and the
 * server - which owns the model - updates the diagram.
 */
(function () {
  'use strict';

  function container(outId) {
    return document.getElementById('graph' + outId);
  }

  window.bscStatus = function (outId, msg) {
    var s = document.getElementById(outId + '_status');
    if (!s) return;
    s.textContent = msg || '';
    s.classList.toggle('bsc-status-active', !!msg);
  };

  window.bscEditorMode = function (outId, mode, type) {
    var el = container(outId);
    if (!el || !el.chart) return;
    var net = el.chart;
    net.disableEditMode();
    if (mode === 'node') {
      el.bscNodeType = type;
      net.addNodeMode();
      window.bscStatus(outId, (type === 'solid' ? 'Solid compartment' : 'Compartment') +
        ': click in the diagram where it should be placed (Esc cancels).');
    } else if (mode === 'edge') {
      el.bscEdgeType = type;
      net.addEdgeMode();
      window.bscStatus(outId, (type === 'perm' ? 'Permeable surface' : 'Fluid flow') +
        ': press on the first compartment and drag to the second one (Esc cancels).');
    } else if (mode === 'delete') {
      var sel = net.getSelection();
      if (sel.nodes.length === 0 && sel.edges.length === 0) {
        window.bscStatus(outId, 'Select a compartment, flow or permeable surface first.');
        return;
      }
      window.bscStatus(outId, '');
      net.deleteSelected();
    } else if (mode === 'fit') {
      net.fit({animation: {duration: 300}});
    }
  };

  // called by the add callbacks once a request has been sent to the server
  window.bscEditorDone = function (outId) {
    window.bscStatus(outId, '');
    var el = container(outId);
    if (el && el.chart) {
      setTimeout(function () { el.chart.disableEditMode(); }, 0);
    }
  };

  // Label of an input with HTML content (e.g. subscripts such as K<sub>media/plastic</sub>).
  // Sent from R with session$sendCustomMessage("bsc-set-label", list(id, html)).
  function registerLabelHandler() {
    if (!window.Shiny || !window.Shiny.addCustomMessageHandler) return false;
    window.Shiny.addCustomMessageHandler('bsc-set-label', function (m) {
      var lab = document.querySelector('label[for="' + m.id + '"]');
      if (lab) lab.innerHTML = m.html;
    });
    return true;
  }
  if (!registerLabelHandler()) {
    document.addEventListener('DOMContentLoaded', registerLabelHandler);
  }

  document.addEventListener('keydown', function (e) {
    if (e.key !== 'Escape') return;
    var editors = document.querySelectorAll('.bsc-editor [id^="graph"]');
    for (var i = 0; i < editors.length; i++) {
      if (editors[i].chart) editors[i].chart.disableEditMode();
    }
    var st = document.querySelectorAll('.bsc-status');
    for (var j = 0; j < st.length; j++) {
      st[j].textContent = '';
      st[j].classList.remove('bsc-status-active');
    }
  });
})();
