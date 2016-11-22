



composable_shiny_panels_app <- function(panels, server, prefix="", title=NULL, top=FALSE) {
    ui <- function(request) div(
        HTML('
<style>
table.dataTable.display tbody td { 
    border-top: 0; 
    border-bottom: 0; 
    border-style: none;
    padding-top: 0; 
    padding-bottom: 0; 
    line-height: 1;
    white-space: nowrap;
}

.selectize-dropdown-content {
    max-height: 50em;
}
</style>
'),
        if (!is.null(title)) titlePanel(title),
        do.call(
            if (top) tabsetPanel else navlistPanel, 
            c(
                list(id=paste0(prefix, "-tabset")),
                if (!top) list(widths=c(2,10),well=FALSE), 
                lapply(panels, call_ui, request)
                )),
        div(style="height: 50em")
    )

    app <- composable_shiny_app(ui, server, title=title)
    app$component_panels <- panels
    app
}

