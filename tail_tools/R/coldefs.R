
#
# Column rendering for DT
#

precision_coldef <- function(target, digits) {
    list(
        targets=target,
        render=htmlwidgets::JS('
        function(data,type,row,meta) {
            var x = parseFloat(data);
            return isNaN(x) ? "" : x.toPrecision(',digits,');
        }')
    )
}

fixed_coldef <- function(target, digits) {
    list(
        targets=target,
        render=htmlwidgets::JS('
        function(data,type,row,meta) {
            var x = parseFloat(data);
            return isNaN(x) ? "" : x.toFixed(',digits,');
        }')
    )
}


r_coldef <- function(target, r_low, r_high) {
    list(
        targets=target,
        className="dt-center",
        render=htmlwidgets::JS('
        function(data,type,row,meta) {
            var r = parseFloat(data);
            var r_low = parseFloat(row[',r_low,']);
            var r_high = parseFloat(row[',r_high,']);
            function tr(x) { return (x+1)*100+1.5; }
            return "<tt style=\\"white-space: pre\\">" + 
                (" "+r.toFixed(2)).slice(-5) + 
                "<span style=\\"font-size: 66%;\\"> [" + (" "+r_low.toFixed(2)).slice(-5) + 
                "," + (" "+r_high.toFixed(2)).slice(-5) + "]</span>" +
                "</tt> " +
                "<svg width=203 height=21 style=\\"vertical-align: middle\\">" +
                "<line x1="+tr(-1)+" y1=0 x2="+tr(-1)+" y2=21 stroke=#000 stroke-width=1 />" +
                "<line x1="+tr( 0)+" y1=0 x2="+tr( 0)+" y2=21 stroke=#bbb stroke-width=1 />" +
                "<line x1="+tr( 1)+" y1=0 x2="+tr( 1)+" y2=21 stroke=#000 stroke-width=1 />" +
                "<line x1="+tr(r_low)+" y1=10 x2="+tr(r_high)+" y2=10 stroke=#000 stroke-width=2 stroke-linecap=butt />" +
                "<circle cx="+tr(r)+" cy=10 r=4 fill=#000 />" +
                "</svg>"
        }')
    )
}

n_coldef <- function(targets, maximum, digits=0, color="#88f") {
    list(
        targets=targets,
        render=htmlwidgets::JS(paste0('
        function(data,type,row,meta) {
            var max = ',maximum,';
            var digits = ',digits,';
            var val = parseFloat(data);
            var size = 21*Math.sqrt((isNaN(val)?0:val)/max);
            var pos = (21-size)/2;
            return "<tt>" + (isNaN(val) ? "" : val.toFixed(digits)) + "</tt> " +
                "<svg width=21 height=21 style=\\"vertical-align: middle\\">" +
                "<rect x="+pos+" y="+pos+" width="+size+" height="+size+" style=\\"fill: ',color,';\\" />"+
                "</svg>";
        }'))
    )
}


bar_coldef <- function(targets, maximum, digits=0, color="#88f") {
    list(
        targets=targets,
        render=htmlwidgets::JS(paste0('
        function(data,type,row,meta) {
            var max = ',maximum,';
            var digits = ',digits,';
            var val = parseFloat(data);
            var size = 42*((isNaN(val)?0:val)/max);
            return "<tt>" + (isNaN(val) ? "" : val.toFixed(digits)) + "</tt> " +
                "<svg width=42 height=21 style=\\"vertical-align: middle\\">" +
                "<rect x=0 y=0 width="+size+" height=21 style=\\"fill: ',color,';\\" />"+
                "</svg>";
        }'))
    )
}

rank_coldef <- function(target, fdr) {
    if (is.null(fdr))
        return(list(targets=target))

    list(
        targets=target,
        render=htmlwidgets::JS('function(data,type,row,meta) {
            var fdr = parseFloat(row[',fdr,']);
            var stars = "";
            var fdr_text, color;
            if (fdr <= 0.01)
                fdr_text = fdr.toExponential(0);
            else
                fdr_text = fdr.toPrecision(1);
            if (fdr <= 0.01)
                color = "#0a0"
            else if (fdr <= 0.05)
                color = "#cc0"
            else 
                color = "#888"
            return "<span style=\\"color:" + color + "\\">" + data + " " +
                "<div title=FDR style=\\"display: inline-block; font-size: 80%; width: 4em;\\">(" + 
                fdr_text + ")</div></span>";
        }')
    )
}


confect_coldef <- function(target, confect, low, high) {
    list(
        targets=target,
        render=htmlwidgets::JS('function(data,type,row,meta) {
            var low=',low,';
            var high=',high,';
            var effect=parseFloat(data);
            var color="#000";
            var from, to, bound;
            if (row[',confect,'] == null) {
                from = low;
                to = high
                color = "#d00";
                bound = "";
            } else {
                var confect = parseFloat(row[',confect,']);
                if (effect >= 0.0) {
                    from = confect;
                    to = high;
                } else {
                    from = low;
                    to = confect;
                }
                bound = " (" + (effect>=0?"&gt;":"&lt;") + (confect>=0?" ":"") + confect.toFixed(2) + ")"
            }
            function tr(x) { return (x-low)/(high-low) * 200 + 1.5; }
            return "<tt style=\\"white-space: pre\\">" + 
                   effect.toFixed(2) + 
                   bound +
                   "</tt>" +
                   "<svg width=203 height=21 style=\\"vertical-align: middle\\">" +
                   "<line x1="+tr(low)+" y1=0 x2="+tr(low)+" y2=21 stroke=#bbb stroke-width=1 />" +
                   "<line x1="+tr(high)+" y1=0 x2="+tr(high)+" y2=21 stroke=#bbb stroke-width=1 />" +
                   "<line x1="+tr(0)+" y1=0 x2="+tr(0)+" y2=21 stroke=#bbb stroke-width=1 />" +
                   "<line x1="+tr(from)+" y1=10 x2="+tr(to)+" y2=10 " + 
                         "stroke="+color+" stroke-width=2 stroke-linecap=butt />" +
                   "<circle cx="+tr(effect)+" cy=10 r=4 fill="+color+" />" +
                   "</svg>"
        }'))
}


