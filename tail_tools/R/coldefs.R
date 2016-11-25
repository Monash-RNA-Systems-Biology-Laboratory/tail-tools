
#
# Column rendering for DT
#

precision_coldef <- function(target, digits) {
    list(
        targets=target,
        render=JS('
        function(data,type,row,meta) {
            var x = parseFloat(data);
            return isNaN(x) ? "" : x.toPrecision(',digits,');
        }')
    )
}

fixed_coldef <- function(target, digits) {
    list(
        targets=target,
        render=JS('
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
        render=JS('
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

n_coldef <- function(targets, maximum, digits=0) {
    list(
        targets=targets,
        render=JS('
        function(data,type,row,meta) {
            var max = ',maximum,';
            var digits = ',digits,';
            var val = parseFloat(data);
            var size = 21*Math.sqrt(val/max);
            var pos = (21-size)/2;
            return "<tt>" + val.toFixed(digits) + "</tt> " +
                "<svg width=21 height=21 style=\\"vertical-align: middle\\">" +
                "<rect x="+pos+" y="+pos+" width="+size+" height="+size+" style=\\"fill: #88f;\\" />"+
                "</svg>";
        }')
    )
}

rank_coldef <- function(target, fdr) {
    if (is.null(fdr))
        return(list(targets=target))

    list(
        targets=target,
        render=JS('function(data,type,row,meta) {
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

