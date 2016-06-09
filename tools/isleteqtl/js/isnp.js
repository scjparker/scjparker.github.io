//
// Copyright 2015 Stephen Parker
//
// Licensed under Version 3 of the GPL or any later version
//

/* global d3, hg19References, hg19ReferenceByName, hg19ReferenceByPosition, chromatinColors, data, footprints, biotypes, eQTLChromatinStateEnrichment, eQTLChromatinStateEnrichmentScaled */

'use strict';

var margin = {top: 10, right: 50, bottom: 100, left: 70};
var plot;
var width;
var height;

var xValue;
var xMap;
var xRange;
var xScale;
var xAxis;
var yValue;
var yMap;
var yRange;
var yScale;
var yAxis;

var zoom;
var tip;
var svg;
var objects;

var biotypeLegendContainer;
var chromatinStateLegendContainer;
var footprintLegendContainer;

var filters = {
    click: {
        chromatinStateFilters: [],
        footprintFilters: [],
        biotypeFilters: []
    },
    hover: {
        chromatinStateFilters: [],
        footprintFilters: [],
        biotypeFilters: []
    }
};

var geneQuery = null;

var consoleEnabled = true;
var snpDetails = {};

var gburlbase = 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=hensley&hgS_otherUserSessionName=20160603_islet_eQTL&position=';

function log(msg) {
    if (consoleEnabled) {
        console.log(msg);
    }
}

function report(evt) {
    var addr;
    var link = evt.currentTarget;
    addr = link.innerText.replace('at', '@').split(/\s/);
    link.href = ['m', 'a', 'i', 'l', 't', 'o', ':'].concat(addr).concat('?subject=Islet%20eQTL%20explorer').join('');
}

function makeReporters() {
    var reporter;
    for (reporter of querySelectorAll('.reporter')) {
        reporter.addEventListener('click', report);
    }
}

var format_number_for_locale = new Intl.NumberFormat(window.navigator.languages || [window.navigator.language || window.navigator.userLanguage], { localeMatcher: 'best fit', maximumFractionDigits: 4 }).format;

d3.selection.prototype.moveToFront = function() {
    return this.each(function() {
        if (this.nextSibling) this.parentNode.appendChild(this);
    });
};

function fillColor(d) {
    var cs = typeof(d) === 'string' ? d : d.cs;
    return chromatinColors[cs].fill;
}

function strokeColor(d) {
    var cs = typeof(d) === 'string' ? d : d.cs;
    return chromatinColors[cs].stroke;
}

function getShape(o) {
    var fp = (typeof(o) === 'object') ? o.fp : o;
    var shape = fp === 'yes' ? 'square' : 'circle';
    return shape;
}

var drawing = false;

function querySelectorAll(selectors, elementOrId) {
    var container;
    selectors = typeof(selectors) == 'string' ? [selectors] : selectors;
    if (elementOrId) {
        container = typeof(elementOrId) == 'string' ? document.getElementById(elementOrId) : elementOrId;
    } else {
        container = document;
    }
    var result = selectors.map(function(selector) {return Array.from(container.querySelectorAll(selector));}).reduce(function(p, c, i, a) {return p.concat(c);}, []);
    return result;
}

function changeClassList(selector, classesToAdd, classesToRemove) {
    querySelectorAll(selector).map(function(el) {
        var classList = el.classList;
        if (classesToAdd.length > 0) {
            classList.add.apply(classList, classesToAdd);
        }
        if (classesToRemove.length > 0) {
            classList.remove.apply(classList, classesToRemove);
        }
    });
}

function getWidth(el) {
    var width = el.offsetWidth;
    var style = window.getComputedStyle(el);
    width = width - parseInt(style.paddingLeft) - parseInt(style.paddingRight);
    return width;
}

function getHeight(el) {
    var height = el.offsetHeight;
    var style = window.getComputedStyle(el);
    height = height - parseInt(style.paddingTop) - parseInt(style.paddingBottom);
    return height;
}

function clickFiltersExist() {
    return (filters.click.biotypeFilters.length > 0 || filters.click.footprintFilters.length > 0 || filters.click.chromatinStateFilters.length > 0);
}

function hoverFiltersExist() {
    return (filters.hover.biotypeFilters.length > 0 || filters.hover.footprintFilters.length > 0 || filters.hover.chromatinStateFilters.length > 0);
}

function filtersExist() {
    return (geneQuery !== null || clickFiltersExist() || hoverFiltersExist());
}

function matchesFilters(filters, filter) {
    return filters.length > 0 && filters.indexOf(filter) !== -1;
}

function classForDot(d) {
    var classes = [
        'dot'
    ];
    classes.push(d.cs);
    classes.push(d.bt);
    classes.push(d.fp);
    return classes.join(' ');
}

function combineArrays(arrays) {
    var combinations = [];
    var filtered = arrays.filter(function(a) {return a.length > 0;}).sort(function(a, b) {return b.length - a.length;});

    if (filtered.length == 1) {
        combinations = filtered[0];
    } else if (filtered.length > 1) {
        combinations = filtered.reduce(function(p, c, i, a) {
            return p.map(function(ps) {
                return c.map(function(cs) {
                    if (Array.isArray(ps)) {
                        return ps.map(function(ps2) {
                            if (Array.isArray(ps2)) {
                                return ps2.map(function(ps3) {
                                    return ps3 + cs;
                                });
                            } else {
                                return ps2 + cs;
                            }
                        });
                    } else {
                        return '.dot' + ps + cs;
                    }
                });
            });
        });
    }
    return combinations.join(',');
}

function filterDots() {
    if (filtersExist()) {
        var dotSelectors = [];

        var biotypeFilters;
        var chromatinStateFilters;
        var footprintFilters;
        var geneFilters = [];
        var allFilters;

        changeClassList('.dot', ['invisible'], ['highlight']);

        if (geneQuery && geneQuery !== '') {
            geneFilters = ['[data-g*="' + geneQuery.toUpperCase() + '"]', '[data-gn*="' + geneQuery.toUpperCase() + '"]'];
        }

        biotypeFilters = filters.click.biotypeFilters.concat(filters.hover.biotypeFilters).map(function(f) {return '[data-bt=' + f + ']';});
        chromatinStateFilters = filters.click.chromatinStateFilters.concat(filters.hover.chromatinStateFilters).map(function(f) {return '[data-cs=' + f + ']';});
        footprintFilters = filters.click.footprintFilters.concat(filters.hover.footprintFilters).map(function(f) {return '[data-fp=' + f + ']';});
        allFilters = [biotypeFilters, chromatinStateFilters, footprintFilters, geneFilters];
        dotSelectors = combineArrays(allFilters);

        if (dotSelectors.length > 0) {
            changeClassList(dotSelectors, ['highlight'], ['invisible']);
        }
    } else {
        changeClassList('.highlight, .invisible', [], ['invisible', 'highlight']);
    }
}

function addFilter(filters, selector) {
    if (filters.indexOf(selector) == -1) {
        filters.push(selector);
    }
}

function removeFilter(filters, selector) {
    var index;
    while ((index = filters.indexOf(selector)) != -1){
        filters.splice(index, 1);
    }
}

function changeFilters(filters, filteringElement) {
    var filter = filteringElement.dataset ? filteringElement.dataset.filter : null;
    if (filter) {
        if (filters.indexOf(filter) == -1) {
            filteringElement.classList.add('selected');
            addFilter(filters, filter);
        } else {
            filteringElement.classList.remove('selected');
            removeFilter(filters, filter);
        }
        return true;
    }
    return false;
}

function makeDotHighlighter(filters) {
    return function(evt) {
        evt.stopPropagation();
        var filter = evt.currentTarget.dataset.filter;
        addFilter(filters, filter);
        requestAnimationFrame(filterDots);
    };
}

function makeDotHighlightResetter(filters) {
    return function(evt) {
        evt.stopPropagation();
        var filter = evt.currentTarget.dataset.filter;
        removeFilter(filters, filter);
        requestAnimationFrame(filterDots);
    };
}

function handleBiotypeClick(evt) {
    evt.stopPropagation();
    if (changeFilters(filters.click.biotypeFilters, evt.currentTarget)) {
        requestAnimationFrame(filterDots);
    }
}

function handleChromatinStateClick(evt) {
    evt.stopPropagation();
    if (changeFilters(filters.click.chromatinStateFilters, evt.currentTarget)) {
        requestAnimationFrame(filterDots);
    }
}

function handleFootprintClick(evt) {
    evt.stopPropagation();
    if (changeFilters(filters.click.footprintFilters, evt.currentTarget)) {
        requestAnimationFrame(filterDots);
    }
}


function setSearchFilter(s) {
    var changed = false;
    if (!s || s === '' && geneQuery) {
        changed = true;
        geneQuery = null;
    } else if (geneQuery === null || s != geneQuery.source) {
        changed = true;
        geneQuery = s; // new RegExp(s, 'ig');
    }
    if (changed) {
        requestAnimationFrame(filterDots);
    }
}

function handleSearchInput(evt) {
    evt.stopPropagation();
    window.setTimeout(function() {setSearchFilter(evt.target.value);}, 100);
}

function reset() {
    filters.click.biotypeFilters = [];
    filters.click.chromatinStateFilters = [];
    filters.click.footprintFilters = [];
    geneQuery = null;
    document.getElementById('search').value = '';
    draw();
}

function clearInput(evt) {
    evt.preventDefault();
    var input = evt.currentTarget.previousElementSibling;
    input.value = '';
    var event = document.createEvent('HTMLEvents');
    event.initEvent('change', true, false);
    input.dispatchEvent(event);
}

function delabel(s) {
    return (typeof(s) === 'string') ? s.replace(/_/g, ' ') : s;
}

function makeEnrichmentBar(cs) {
    var change = eQTLChromatinStateEnrichment[cs];
    var scaledChange = eQTLChromatinStateEnrichmentScaled[cs];
    var x = 4 * (scaledChange < 0 ? (1 + scaledChange) : 1);
    var width = Math.abs(scaledChange) * 4;
    var height = '1.25em';
    var y = '0.1em';
    var textX = change < 0 ? '0.5em' : '5.5em';
    var content = '<div class="enrichment">' +
        '<svg>' +
        '<rect y="' + y + '" height="' + height + '" x="' + x + 'em" width="' + width + 'em" fill="#ddd" stroke="#999" stroke-width="0.2px"/>' +
        '<line x1="4em" x2="4em" y1="0.1em" y2="1.35em" stroke="#777" stroke-width="1px"/>' +
        '<rect x="0.1em" y="0.1em" height="1.25em" width="8em" stroke-width="0.5px" stroke="#000" fill="transparent"/>' +
        '<text x="' + textX + '" y="1.05em">' + change + '</text>' +
        '<text x="3.75em" y="1.05em">0</text>' +
        '</svg></div>';
    return content;
}

function makeChromatinStateContent(cs) {
    return '<span class="legendColor"><svg><rect x="0.1em" y="0.1em" width="1em" height="1em" fill="' + fillColor(cs) + '" stroke="' + strokeColor(cs) + '" stroke-width="0.5px"/></svg></span><span class="label">' + delabel(cs) + '</span>';
}

var footPrintContent = {
    'yes': '<span class="legendColor"><svg class="rotate45"><rect x="0.1em" y="0.1em" width="0.85em" height="0.85em" fill="#999"/></svg></span>yes',
    'no': '<span class="legendColor"><svg><circle cx="0.5em" cy="0.5em" r="0.5em" fill="#999"/></svg></span>no'
};

function makeFootprintContent(o) {
    var fp = (typeof(o) === 'object') ? o.fp : o;
    return footPrintContent[fp];
}

function placeDot(d) {
    var tx = 'translate(' + xMap(d) + ',' + yMap(d) + ')';
    if (d.fp === 'yes') {
        tx += ' rotate(45)';
    }
    return tx;
}

function makeSNPSortKey(snpID) {
    var chr = Number(snpID.split(':')[0]);
    var pos = Number(snpID.split(':')[1].split('_')[0]);
    var alleles = snpID.split('_')[1];
    return [chr, pos, alleles];
}

function sortSNPs(a, b) {
    var i;
    var cmp;
    var k1 = makeSNPSortKey(a);
    var k2 = makeSNPSortKey(b);
    for (i = 0; i < k1.length; i++) {
        cmp = k1[i] - k2[i];
        if (cmp !== 0) {
            break;
        }
    }
    return cmp;
}

function removeSNPDetails(evt) {
    evt.stopPropagation();
    var tr = evt.currentTarget.parentNode.parentNode;
    delete snpDetails[tr.firstChild.innerHTML];
    requestAnimationFrame(makeSNPDetailsTable);
}

function makeSNPDetailsTable() {
    var tr;
    var td;
    var remover;

    var tbody = document.getElementById('snpDetailsBody');

    var snps = Object.keys(snpDetails).sort(sortSNPs);
    if (snps.length > 0) {
        tbody.innerHTML = '';
        for (var snp of snps) {
            tr = document.createElement('tr');
            for (var field of snpDetails[snp]) {
                td = document.createElement('td');
                td.innerHTML = field;
                tr.appendChild(td);
            }
            td = document.createElement('td');
            remover = document.createElement('button');
            remover.innerHTML = 'Remove';
            remover.addEventListener('click', removeSNPDetails, true);
            td.appendChild(remover);
            tr.appendChild(td);
            tbody.appendChild(tr);
        }
    } else {
        tbody.innerHTML = '<td colspan="9"><span class="instruction">Click on a SNP in the plot to show its details here.</span></td>';
    }
}

function addSNPDetails(d, i) {
    d3.event.preventDefault();
    hideTooltip(d, i);
    var gbwindow = d.chr + ':' + (d.pos - 50000) + '-' + (d.pos + 49999);
    var fields = [
        d.snp,
        '<a target="_blank" href="http://www.ensembl.org/id/' + d.g + '">' + d.g + '</a>',
        d.gn,
        '<a target="_blank" href="' + gburlbase + gbwindow + '">' + d.chr + ':' + format_number_for_locale(d.pos) + '</a>',
        delabel(d.bt),
        makeChromatinStateContent(d.cs),
        makeFootprintContent(d),
        d.lp,
        d.a1 + '/' + d.a2
    ];
    snpDetails[d.snp] = fields;
    if (Object.keys(snpDetails).length > 0) {
        requestAnimationFrame(makeSNPDetailsTable);
    }
}

function showTooltip(d, i) {
    d3.select(d3.event.target).moveToFront();
    d3.event.target.classList.add('opaque');
    tip.show(d, i);
}

function hideTooltip(d, i) {
    d3.event.target.classList.remove('opaque');
    tip.hide(d, i);
}

function draw() {
    if (drawing) {
        return;
    }
    drawing = true;

    var li;
    plot = document.getElementById('plot');
    plot.innerHTML = '';
    d3.select('.d3-tip').remove();
    document.getElementById('chromatinStateLegend').innerHTML = '';
    document.getElementById('footprintLegend').innerHTML = '';
    document.getElementById('biotypeLegend').innerHTML = '';

    width = Math.max(getWidth(plot) - margin.left - margin.right);
    var controlHeight = getHeight(document.getElementById('control'));
    height = Math.max(getWidth(plot) * 0.4, controlHeight);

    xValue = function(d) { return d.pos + hg19ReferenceByName[d.chr].offset;};
    xMap = function(d) { return xScale(xValue(d));};
    xRange = [-50000000, hg19References[hg19References.length - 1].offset + hg19References[hg19References.length - 1].length + 50000000];

    xScale = d3.scale.linear()
        .domain(xRange)
        .range([0, width]);

    xAxis = d3.svg.axis()
        .orient('bottom')
        .scale(xScale)
        .tickValues(d3.set(hg19References.map(function(r) {return r.offset;})).values())
        .tickFormat(function(r) {return hg19ReferenceByPosition[r].name;});

    yValue = function(d) { return d.lp;};
    yMap = function(d) { return yScale(yValue(d));};
    yRange = d3.extent([-1, d3.max(data, function(d) { return d.lp; }) + 5]);

    yScale = d3.scale.linear().domain(yRange).range([height, 0]);

    yAxis = d3.svg.axis().orient('left').scale(yScale);

    zoom = d3.behavior.zoom()
        .scaleExtent([1, 32])
        .x(xScale)
        .y(yScale)
        .on('zoom', function() {
            d3.select('.x.axis').call(xAxis);
            d3.select('.y.axis').call(yAxis);
            svg.selectAll('.dot').attr('transform', placeDot);
        });

    tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([-10, 10])
        .direction('se')
        .html(
            function(d) {
                return '<div><h3>SNP ' + d.snp + '</h3>' +
                    '<table><tbody>' +
                    '<tr><th>Gene</th><td>' + d.g + '</td></tr>' +
                    '<tr><th>Gene name</th><td>' + d.gn + '</td></tr>' +
                    '<tr><th>Position</th><td>' + d.chr + ':' + d.pos + '</td></tr>' +
                    '<tr><th>Biotype</th><td>' + delabel(d.bt) + '</td></tr>' +
                    '<tr><th>Chromatin state</th><td>' + makeChromatinStateContent(d.cs) + '</td></tr>' +
                    '<tr><th>ATAC-seq footprint?</th><td>' + makeFootprintContent(d) + '</td></tr>' +
                    '<tr><th>eQTL -log<sub>10</sub>(p-value)</th><td>' + d.lp + '</td></tr>' +
                    '<tr><th>Allele 1/2</th><td>' + d.a1 + '/' + d.a2 + '</td></tr>' +
                    '<tbody></table></div>';
            }
        );

    document.getElementById('plot').innerHTML = '';

    var totalWidth = width + margin.left + margin.right;
    var totalHeight = height + margin.top + margin.bottom;

    svg = d3.select('#plot').append('svg')
        .attr('width', totalWidth)
        .attr('height', totalHeight)
        .attr('viewBox', '0 0 ' + totalWidth + ' ' + totalHeight)
        .attr('preserveAspectRatio', 'xMinYMin')
        .append('g')
        .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')
        .call(zoom);

    svg.call(tip);

    svg.append('rect')
        .attr('fill', 'rgba(255, 255, 255, 0)')
        .attr('width', width)
        .attr('height', height);

    svg.append('g')
        .attr('class', 'x axis')
        .attr('transform', 'translate(0,' + height + ')')
        .call(xAxis)
        .append('text')
        .attr('class', 'label')
        .attr('x', width * 0.5)
        .attr('y', 40)
        .style('text-anchor', 'middle')
        .text('Genomic Location');

    svg.append('g')
        .attr('class', 'y axis')
        .call(yAxis)
        .append('text')
        .attr('class', 'label')
        .attr('transform', 'rotate(-90), translate(0, -20)')
        .attr('x', -1 * height / 2)
        .attr('y', -35)
        .style('text-anchor', 'middle')
        .text('eQTL -log')
        .append('tspan')
        .attr('class', 'label')
        .attr('transform', 'rotate(-90), translate(0, -20)')
        .attr('dy', '.70em')
        .style('font-size', '0.8em')
        .text('10')
        .append('tspan')
        .attr('class', 'label')
        .attr('transform', 'rotate(-90), translate(0, -20)')
        .attr('dy', '-0.70em')
        .style('font-size', '1.2em')
        .style('text-anchor', 'middle')
        .text('(p-value)');

    objects = svg.append('svg')
        .classed('objects', true)
        .attr('width', width)
        .attr('height', height);

    objects.selectAll('.dot')
        .data(data)
        .enter().append('path')
        .attr('d', d3.svg.symbol().type(getShape))
        .attr('class', classForDot)
        .attr('id', function(d) {return d.id;})
        .attr('data-g', function(d) {return d.g;})
        .attr('data-gn', function(d) {return d.gn;})
        .attr('data-cs', function(d) {return d.cs;})
        .attr('data-bt', function(d) {return d.bt;})
        .attr('data-fp', function(d) {return d.fp;})
        .attr('transform', placeDot)
        .style('fill', fillColor)
        .style('stroke', strokeColor)
        .on('mouseover', showTooltip)
        .on('mouseout', hideTooltip)
        .on('touchstart', showTooltip)
        .on('touchmove', hideTooltip)
        .on('click', addSNPDetails)
        .on('touchend', addSNPDetails);

    chromatinStateLegendContainer = document.getElementById('chromatinStateLegend');
    li = document.createElement('li');
    li.classList.add('header');
    li.innerHTML = '<span class="legendColor">State</span><span class="enrichment">eQTL log2 fold enrichment</span>';
    chromatinStateLegendContainer.appendChild(li);

    Object.keys(chromatinColors).map(function(chromatinState) {
        li = document.createElement('li');
        if (matchesFilters(filters.click.chromatinStateFilters, chromatinState)) {
            li.classList.add('selected');
        }
        li.dataset.filter = chromatinState;
        li.innerHTML = makeChromatinStateContent(chromatinState) + makeEnrichmentBar(chromatinState);
        li.addEventListener('mouseover', makeDotHighlighter(filters.hover.chromatinStateFilters), true);
        li.addEventListener('mouseleave', makeDotHighlightResetter(filters.hover.chromatinStateFilters), true);
        li.addEventListener('click', handleChromatinStateClick, true);
        chromatinStateLegendContainer.appendChild(li);
    });

    footprintLegendContainer = document.getElementById('footprintLegend');

    footprints.map(function(footprint) {
        li = document.createElement('li');
        if (matchesFilters(filters.click.footprintFilters, footprint)) {
            li.classList.add('selected');
        }
        li.dataset.filter = footprint;
        li.innerHTML = makeFootprintContent(footprint);
        li.addEventListener('mouseover', makeDotHighlighter(filters.hover.footprintFilters), true);
        li.addEventListener('mouseleave', makeDotHighlightResetter(filters.hover.footprintFilters), true);
        li.addEventListener('click', handleFootprintClick, true);
        footprintLegendContainer.appendChild(li);
    });

    biotypeLegendContainer = document.getElementById('biotypeLegend');

    biotypes.map(function(biotype) {
        li = document.createElement('li');
        if (matchesFilters(filters.click.biotypeFilters, biotype)) {
            li.classList.add('selected');
        }
        li.dataset.filter = biotype;
        li.innerHTML = delabel(biotype);
        li.addEventListener('mouseover', makeDotHighlighter(filters.hover.biotypeFilters), true);
        li.addEventListener('mouseleave', makeDotHighlightResetter(filters.hover.biotypeFilters), true);
        li.addEventListener('click', handleBiotypeClick, true);
        biotypeLegendContainer.appendChild(li);
    });

    requestAnimationFrame(filterDots);
    drawing = false;
}

window.addEventListener('resize', function() {requestAnimationFrame(draw);});
window.addEventListener('orientationchange', function() {requestAnimationFrame(draw);});
document.getElementById('search').addEventListener('change', handleSearchInput);
document.getElementById('search').addEventListener('keyup', handleSearchInput);
document.getElementById('reset').addEventListener('click', reset);
for (var clearButton of querySelectorAll('.inputbox .cleaner')) {
    clearButton.addEventListener('click', clearInput, true);
}
document.getElementById('search').value = '';
makeReporters();
draw();
