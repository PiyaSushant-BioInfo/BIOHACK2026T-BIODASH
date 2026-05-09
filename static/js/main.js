/* TBAnalytica — Frontend Logic */

// ── Section visibility ────────────────────────────────────────────────────

function showSection(id) {
    ['home-section', 'form-section', 'progress-section', 'error-section'].forEach(function(s) {
        var el = document.getElementById(s);
        if (el) el.classList.add('hidden');
    });
    var target = document.getElementById(id);
    if (target) target.classList.remove('hidden');
}

function showHome() {
    showSection('home-section');
    hideAllForms();
}

function hideAllForms() {
    var forms = document.querySelectorAll('.form-card');
    forms.forEach(function(f) { f.classList.add('hidden'); });
}

function showForm(mode) {
    showSection('form-section');
    hideAllForms();
    var formEl = document.getElementById('form-' + mode);
    if (formEl) formEl.classList.remove('hidden');
}

// ── New variant: toggle accession / sequence ──────────────────────────────

function toggleNewInput(type) {
    var accGroup = document.getElementById('accession-group');
    var seqGroup = document.getElementById('sequence-group');
    if (type === 'accession') {
        accGroup.classList.remove('hidden');
        seqGroup.classList.add('hidden');
    } else {
        accGroup.classList.add('hidden');
        seqGroup.classList.remove('hidden');
    }
}

// ── Form submissions ─────────────────────────────────────────────────────

function submitKnown(e) {
    e.preventDefault();
    var variantId = document.getElementById('known-variant').value;
    var patientId = document.getElementById('known-patient').value;
    if (!variantId) { alert('Please select a lineage or variant.'); return false; }
    runJob({
        mode: 'known',
        variant_id: variantId,
        patient_id: patientId
    });
    return false;
}

function submitNew(e) {
    e.preventDefault();
    var inputType = document.querySelector('input[name="input_type"]:checked').value;
    var accession = document.getElementById('new-accession').value.trim();
    var sequence = document.getElementById('new-sequence').value.trim();
    var geneName = document.getElementById('new-gene').value;
    var patientId = document.getElementById('new-patient').value.trim();

    if (inputType === 'accession' && !accession) {
        alert('Please enter an NCBI accession number.');
        return false;
    }
    if (inputType === 'sequence' && !sequence) {
        alert('Please enter a sequence.');
        return false;
    }

    runJob({
        mode: 'new',
        accession: inputType === 'accession' ? accession : '',
        sequence: inputType === 'sequence' ? sequence : '',
        gene_name: geneName,
        patient_id: patientId
    });
    return false;
}

function submitCompare(e) {
    e.preventDefault();
    var first = document.getElementById('cmp-first').value.trim();
    var second = document.getElementById('cmp-second').value.trim();
    if (!first || !second) { alert('Please enter both variant IDs.'); return false; }
    runJob({
        mode: 'compare',
        variant_id: first,
        second_input: second
    });
    return false;
}

function startDemo() {
    runJob({ mode: 'demo' });
}

// ── Job execution + polling ──────────────────────────────────────────────

var currentJobId = null;
var pollTimer = null;

function runJob(params) {
    showSection('progress-section');

    fetch('/run-analysis', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(params)
    })
    .then(function(res) { return res.json(); })
    .then(function(data) {
        if (data.error) {
            showError(data.error);
            return;
        }
        currentJobId = data.job_id;
        renderProgressSteps(data.steps, []);
        startPolling(data.job_id);
    })
    .catch(function(err) {
        showError('Failed to start analysis: ' + err.message);
    });
}

function renderProgressSteps(steps, completed) {
    var container = document.getElementById('progress-steps');
    var html = '';
    var activeFound = false;

    for (var i = 0; i < steps.length; i++) {
        var step = steps[i];
        var isDone = completed.indexOf(step) !== -1;
        var isActive = false;

        if (!isDone && !activeFound) {
            isActive = true;
            activeFound = true;
        }

        var cls = isDone ? 'done' : (isActive ? 'active' : 'pending');
        var iconCls = isDone ? 'step-icon-done' : (isActive ? 'step-icon-active' : 'step-icon-pending');
        var icon = isDone ? '&#10003;' : (isActive ? '&#10227;' : '&#9675;');

        html += '<div class="progress-step ' + cls + '">';
        html += '<span class="step-icon ' + iconCls + '">' + icon + '</span>';
        html += '<span>' + step + '</span>';
        html += '</div>';
    }

    container.innerHTML = html;

    var spinner = document.getElementById('spinner');
    var allDone = completed.length >= steps.length;
    spinner.style.display = allDone ? 'none' : 'block';
}

function startPolling(jobId) {
    if (pollTimer) clearInterval(pollTimer);

    pollTimer = setInterval(function() {
        fetch('/status/' + jobId)
        .then(function(res) { return res.json(); })
        .then(function(data) {
            if (data.error && data.status !== 'running') {
                clearInterval(pollTimer);
                showError(data.error);
                return;
            }

            renderProgressSteps(data.steps, data.progress);

            // Show redirect notification if the server switched pipelines
            if (data.notification) {
                var noteEl = document.getElementById('progress-notification');
                if (noteEl) {
                    noteEl.textContent = data.notification;
                    noteEl.classList.remove('hidden');
                }
            }

            if (data.status === 'complete') {
                clearInterval(pollTimer);
                setTimeout(function() {
                    window.location.href = '/results/' + jobId;
                }, 600);
            } else if (data.status === 'error') {
                clearInterval(pollTimer);
                showError(data.error || 'Unknown error occurred.');
            }
        })
        .catch(function(err) {
            clearInterval(pollTimer);
            showError('Lost connection: ' + err.message);
        });
    }, 800);
}

function showError(msg) {
    showSection('error-section');
    document.getElementById('error-message').textContent = msg;
}

// ── Results page: tabs ───────────────────────────────────────────────────

document.addEventListener('DOMContentLoaded', function() {
    var tabs = document.querySelectorAll('.tab');
    tabs.forEach(function(tab) {
        tab.addEventListener('click', function() {
            var target = this.getAttribute('data-tab');

            document.querySelectorAll('.tab').forEach(function(t) {
                t.classList.remove('active');
            });
            this.classList.add('active');

            document.querySelectorAll('.tab-panel').forEach(function(p) {
                p.classList.remove('active');
            });
            var panel = document.getElementById('tab-' + target);
            if (panel) panel.classList.add('active');
        });
    });
});
