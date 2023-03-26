#' @title  scrnaVis for scRNA-seq visualized elegantly
#' @details  This package provides multiple ways to interactively display your single-cell data
#' @param object a seurat object which has run seurat pipeline
#' @param markers  a gene list of cell-specific expression
#' @import irGSEA gginnards msigdbr egg ggsci ggplot2 ComplexHeatmap dplyr Seurat
#' @import ggstatsplot shinydashboard shiny Cairo
#' @return shiny
#' @export
#' @examples
#' \donttest{
#'   pbmc=readRDS(system.file("data","pbmc.rda",package="scrnaVis"))
#'   scrnaVis(object=pbmc,markers=c("CD3E","CD3G","CD3D"))
#' }
scrnaVis <- function(object=NULL, markers=NULL) {
  gene_markers <- unique(markers)
  ident <- colnames(object@meta.data)
  getsetdb <- c("H", "C2", "C5", "C8")
  geneset <- c(
    msigdbr::msigdbr("Homo sapiens", category = "H") %>% dplyr::distinct(gs_name) %>% .$gs_name,
    msigdbr::msigdbr("Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% dplyr::distinct(gs_name) %>% .$gs_name,
    msigdbr::msigdbr("Homo sapiens", category = "C5", subcategory = "BP") %>% dplyr::distinct(gs_name) %>% .$gs_name,
    msigdbr::msigdbr("Homo sapiens", category = "C8") %>% dplyr::distinct(gs_name) %>% .$gs_name
  )
  gene <- row.names(object)

# adapt to height
  heigth <- function(x) {
    if(length(x) %in% c(0:8)) {
      return("400px")
    }else if(length(x) %in% c(9:16)) {
      return("600px")
    }else if(length(x) > 16) {
      return("800px")
    }
  }
  ############################## > Header  ##################################
    header <- dashboardHeader(title = "ScRNA-seq Vis")

  ############################## > sidebar  ##################################
    sidebar <- dashboardSidebar(
      sidebarMenu(
        menuItem(
          tabName = "markerVis",
          text = "Marker Visualization",
          icon = icon("star")
        ),
        menuItem(
          tabName = "enrichVis",
          text = "Gene Set Exploration",
          icon = icon("bolt")
        )
      )
    )

  ############################## > body  ##################################


    body = dashboardBody(
      tabItems(
        tabItem(
          tabName = "markerVis",
        ############################## > submit parameter  ##################################
          box(
            width = 12,
            div(
              style = "display: inline-block;vertical-align:top; width: 49%;",
              selectizeInput(
                "select_ident","Ident",
                choice = ident,selected = NULL,multiple = F,
                options = list(
                  delimiter = ',',
                  create = I(
                    "function(input, callback){
                        return {
                        value: input,
                        text: input
                        };
                     }"
                  )
                )
              )
            ),
            br(),
            div(
              style = "display: inline-block;vertical-align:top; width: 49%;",
              actionButton(
                "start",
                HTML("<b>Start Analysis</b>"),
                style = "color:red",
                icon = icon("paper-plane")
              )
            )
        ),
        ############################## > TSNE  ##################################
        fluidRow(
          column(6,
            box(
              width = 12,height = "800px",
              div(
                style = "display: inline-block;vertical-align:top; width: 49%;",
                downloadButton(
                  "download_TSNEPlot",
                  HTML("<b>Download</b>"),
                  icon = icon("download")
                )
              ),
              br(),
              plotOutput("TSNE_Fig", width = "80%", height = "400px")
            )
          ),
        ############################## > Statistical proportion  ##################################
          column(6,
            box(
              width = 12,height = "800px",
              div(
                style = "display: inline-block;vertical-align:top; width: 49%;",
                selectizeInput(
                  "select_groupBy",
                  "GroupBy",
                  choice = ident,selected = "orig.ident",multiple = F,
                  options = list(
                    delimiter = ',',
                    create = I(
                      "function(input, callback) {
                         return {
                            value: input,
                            text: input
                                };
                       }"
                    )
                  )
                )
              ),
              br(),
              div(
                style = "display: inline-block;vertical-align:top; width: 49%;",
                downloadButton(
                  "download_PropPlot",
                  HTML("<b>Download</b>"),
                  icon = icon("download")
              )
            ),
            br(),
            plotOutput("PropPlot", height = "400px")
          )
         )
        ),
        fluidRow(
          column(6,
            tabBox(
              width = 12,height = "800px",
              tabPanel(
                HTML("<b>DotPlot</b>"),
                br(),
                div(
                  style = "display: inline-block;vertical-align:top; width: 49%;",
                  downloadButton(
                    "download_DotPlot",
                    HTML("<b>Download</b>"),
                    icon = icon("download")
                  )
                ),
                br(),
                plotOutput("DotPlot", width = "80%",height=heigth(gene_markers))
              ),
              tabPanel(
                HTML("<b>HeatPlot</b>"),
                br(),
                div(
                  style = "display: inline-block;vertical-align:top; width: 49%;",
                  downloadButton(
                    "download_HeatmapPlot",
                    HTML("<b>Download</b>"),
                    icon = icon("download")
                  )
                ),
                br(),
                plotOutput("HeatmapPlot", width = "80%",height=heigth(gene_markers))
              )
            )
          ),
        ############################## > FeaturePlot  ##################################
          column(6,
            box(
              width = 12,height = "800px",
              div(
                style = "display: inline-block;vertical-align:top; width: 100%;",
                selectizeInput(
                  "select_gene","Genes Input",
                  choice = gene,selected = NULL,multiple = T ,
                  options = list(
                    delimiter = ',',
                    create = I(
                      "function(input, callback){
                          return {
                            value: input,
                            text: input
                                  };
                       }"
                    )
                  )
                )
              ),
              br(),
              div(
                style = "display: inline-block;vertical-align:top; width: 19%;",
                actionButton(
                  "submit_gene",
                  HTML("<b>Submit</b>"),
                  icon = icon("paper-plane")
                )
              ),
              div(
                style = "display: inline-block;vertical-align:top; width: 50%;",
                downloadButton(
                  "download_FeaturePlot",
                  HTML("<b>Download</b>"),
                  icon = icon("download")
                )
              ),
              br(),
              plotOutput("FeaturePlot", width = "80%", height ="400px")
            )
          )
        )
      ),
      ############################## > GSVA  ##################################
      tabItem(
        tabName = "enrichVis",
        box(
          width = 12,
          br(),
          div(
            style = "display: inline-block;vertical-align:top; width: 34%;",
            selectInput("select_geneset_db", "Gene Sets Categories",choices = getsetdb)
          ),
          div(
            style = "display: inline-block;vertical-align:top; width: 34%;",
            selectInput("select_geneset_sub_db", "Gene Sets Sub Categories",choices = NULL)
          ),
          br(),
          div(
            style = "display: inline-block;vertical-align:top; width: 19%;",
            actionButton(
              "submit_gene_db",
              HTML("<b>Submit</b>"),
              width = 200,
              icon = icon("paper-plane")
            )
          )
        ),
        fluidRow(
          column(6,
            box(
              width = 12,
              height = "600px",
              HTML("<b>GSVA HeatmapPlot</b>"),
              br(),
              div(
                style = "display: inline-block;vertical-align:top; width: 49%;",
                downloadButton(
                  "download_Enri_HeatmapPlot",
                  HTML("<b>Download</b>"),
                  icon = icon("download")
                )
              ),
              br(),
              plotOutput("Enri_HeatmapPlot", width = "100%", height = "400px")
            )
          ),
          column(6,
            box(
              width = 12,
              div(
                style = "display: inline-block;vertical-align:top; width: 100%;",
                selectizeInput(
                  "select_geneset","Gene Sets",
                  choice = geneset,selected = NULL,multiple = F ,
                  options = list(
                    delimiter = ',',
                    create = I(
                      "function(input, callback){
                         return {
                           value: input,
                           text: input
                         };
                       }"
                    )
                  )
                )
              ),
              div(
                style = "display: inline-block;vertical-align:top; width: 49%;",
                actionButton(
                  "submit_geneset",
                  HTML("<b>Submit</b>"),
                  icon = icon("paper-plane")
                )
              ),
              br(),
              br(),
              tabBox(
                width = 12,
                tabPanel(
                  HTML("<b>VlnPlot</b>"),
                  br(),
                  plotOutput("Enri_VlnPlot", width = "80%", height ="400px")
                ),
                tabPanel(
                  HTML("<b>FeaturePlot</b>"),
                         br(),
                         plotOutput("Enri_FeaturePlot", width = "80%", height ="400px")
                )
              )
            )
          )
        )
      )
    )
  )
  ui <- dashboardPage(header, sidebar, body)

  server <- function(input, output, session) {
    observeEvent(input$start, {
    #显示运行
      showNotification(
        h1("START RUNNING", style = "color:red"),
        duration = 20,
        closeButton = F,
        action = a(href = "javascript:location.reload();")
      )
      ############################## > TSNE  ##################################
      p1 <- TSNEPlot(object,group.by = input$select_ident,label = T) + ggtitle("")

      output$TSNE_Fig <- renderPlot({
        p1
      })

      output$download_TSNEPlot <- downloadHandler(
        filename = function() {
          "TNSEPlot.pdf"
        },
        content = function(file) {
          ggsave(p1, filename = file)
        }
      )
      ############################## > Proportion statistics  ##################################
      PropPlot <- function(group, groupBy) {
        Prop_data <- object@meta.data %>%
          dplyr::select({{group}}, {{groupBy}}) %>%
          dplyr::rename(group = {{group}}) %>%
          dplyr::rename(groupBy = {{groupBy}})

      ## plot proportion
        Prop_fig <- ggbarstats(
          data = Prop_data,
          x = group,
          y = groupBy,
          palette = 'category20c_d3',
          package = "ggsci",
          results.subtitle = F,
          bf.message = F,
          proportion.test = F,
          label.args = list(
            size = 2,
            fill = 'white',
            alpha = 0.85,
            family = 'Arial',
            fontfacet = 'bold'
          ),
          perc.k = 2,
          tilte = '',
          xlab = '',
          legend.title = "Class",
          ggtheme = ggpubr::theme_pubclean()
          ) +
          theme(
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = 'black', lineend = 'round'),
            legend.position = 'right',
            axis.text.x = element_text(size = 15,color = 'black',family = 'Arial'),
            axis.text.y = element_text(size = 15,color = 'black',family = 'Arial'),
            legend.text = element_text(family = 'Arial',size = 10,color = 'black'),
            legend.title = element_text(family = 'Arial',size = 13,color = 'black')
          )

        delete_layers(x = Prop_fig, match_type = 'GeomText')
      }

      p2 <- PropPlot(input$select_ident, input$select_groupBy)

      output$PropPlot <- renderPlot({
        PropPlot(input$select_ident, input$select_groupBy)
      })

      output$download_PropPlot <- downloadHandler(
        filename = function() {
          "PropPlot.pdf"
        },
        content = function(file) {
          ggsave(p2, filename = file)
        }
      )
      ############################## > marker Visualization  ##################################
      # Dotplot
      p3 <- DotPlot(object,features = gene_markers,group.by = input$select_ident) +
        coord_flip() + theme_bw() +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size = 12),
          axis.text.y = element_text(face = "italic", size = 12)
        ) +
        scale_color_gradientn(
          values = seq(0, 1, 0.2),
          colors = c("#330066", "#336699", "#66CC66", "#FFCC33")
        ) +
        labs(x = NULL, y = NULL)

      output$DotPlot <- renderPlot({
        p3
      })
      #download
      output$download_DotPlot = downloadHandler(
        filename = function() {
          "DotPlot.pdf"
        },
        content = function(file) {
          ggsave(p3, filename = file)
        }
      )
      # heatmapPlot
      if (is.null(object@assays$SCT)) {
        heatmap_AveE <- as.data.frame(
          AverageExpression(
            object,features = gene_markers,verbose = F,
            slot = "data",group.by = input$select_ident) %>% .$RNA
        )
        bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 2, by = 0.01))
        heatmap <- ComplexHeatmap::pheatmap(
          heatmap_AveE,
          cluster_cols = F,cluster_rows = F,show_colnames = T,show_rownames = T,border = T,
          #border_color = "white",
          color = c(
            colorRampPalette(colors = c("#2166ac", "#f7fbff"))(length(bk) / 2),
            colorRampPalette(colors = c("#f7fbff", "#b2182b"))(length(bk) / 1)
          ),
          breaks = bk,scale = "row",legend_breaks = seq(-1, 2, 1),name = "Exp"
        )
        heatmap@row_names_param$gp <- grid::gpar(fontface = "italic")
      } else {
          heatmap_AveE <- as.data.frame(
            AverageExpression(
              object,features = gene_markers,
              verbose = F,slot = "data",group.by = input$select_ident) %>% .$SCT
          )
          bk <- c(seq(-1, -0.1, by = 0.01), seq(0, 2, by = 0.01))
          heatmap <- ComplexHeatmap::pheatmap(
            heatmap_AveE,
            cluster_cols = F,cluster_rows = F,show_colnames = T,show_rownames = T,border = T,
              #border_color = "white",
              color = c(
                colorRampPalette(colors = c("#2166ac", "#f7fbff"))(length(bk) / 2),
                colorRampPalette(colors = c("#f7fbff", "#b2182b"))(length(bk) / 1)
              ),
              breaks = bk,scale = "row",legend_breaks = seq(-1, 2, 1),name = "Exp"
            )
          heatmap@row_names_param$gp <- grid::gpar(fontface = "italic")
        }

      output$HeatmapPlot <- renderPlot({
        heatmap
      })
      # download
      output$download_HeatmapPlot <- downloadHandler(
        filename = function() {
          "HeatmapPlot.pdf"
        },
        content = function(file) {
          ggsave(heatmap, filename = file,width=10,height=10)
        }
      )
      # featurePlot
      observeEvent(input$submit_gene, {
        p3 <- FeaturePlot(
              object,features = input$select_gene,
              cols = c('#7F7F7FFF', 'blue', "#D62728FF"),
              pt.size = 0.5,reduction = "tsne"
          ) &
          egg::theme_article() &
          theme(plot.title  = element_text(face = "italic"))

        output$FeaturePlot <- renderPlot({
          p3
        })
        #download
        output$download_FeaturePlot <- downloadHandler(
          filename = function() {
            "FeaturePlot.pdf"
          },
          content = function(file) {
            ggsave(p3, filename = file)
          }
        )
      })
          ############################## > GSVA  ##################################
      observe({
        x <- input$select_geneset_db
        if (x %in% c("H", "C8")) {
          updateSelectInput(session,"select_geneset_sub_db",choices = NULL,selected = NULL)
        } else if (x == "C2") {
            y = "CP:KEGG"
            updateSelectInput(session, "select_geneset_sub_db",choices = y)
        } else if (x == "C5") {
            y = "BP"
            updateSelectInput(session, "select_geneset_sub_db",choices = y)
            }
      })

      observeEvent(input$submit_gene_db, {
        if (is.null(object@assays$SCT)) {
          GSEA <- irGSEA::irGSEA.score(
            object = object,
            assay = "RNA",slot = "data",seeds = 123,ncores = 8,
            min.cells = 3,min.feature = 0,custom = F,geneset = NULL,
            msigdb = T,
            species = "Homo sapiens",
            category = input$select_geneset_db,
            subcategory = input$selcet_geneset_sub_db,
            category = "H",geneid = "symbol",
            method = "AUCell",aucell.MaxRank = NULL,
            ucell.MaxRank = NULL,kcdf = 'Gaussian'
          )
        } else {
            GSEA <- irGSEA::irGSEA.score(
              object = object,
              assay = "SCT",slot = "data",seeds = 123,ncores = 8,
              min.cells = 3,min.feature = 0,custom = F,geneset = NULL,
              msigdb = T,species = "Homo sapiens",
              category = input$select_geneset_db,
              subcategory = input$selcet_geneset_sub_db,
              geneid = "symbol",method = "AUCell",
              aucell.MaxRank = NULL,ucell.MaxRank = NULL,
              kcdf = 'Gaussian'
            )
          }
      result.dge <- irGSEA::irGSEA.integrate(
        object = GSEA,
        group.by = input$select_ident,
        metadata = NULL,
        col.name = NULL,
        method = "AUCell"
      )
            # 热图
      p4 <- irGSEA::irGSEA.heatmap(
        object = result.dge,
        method = "AUCell",
        top = 30,
        show.geneset = NULL
      )
      output$Enri_HeatmapPlot <- renderPlot({
        p4
      })
      output$download_Enri_HeatmapPlot <- downloadHandler(
        filename = function() {
          "Enriment_HeatmapPlot.pdf"
        },
        content = function(file) {
          ggsave(p4, filename = file,width=10,height=10)
        }
      )

      observeEvent(input$submit_geneset, {
        # featurePlot
        output$Enri_FeaturePlot <- renderPlot(
          {
          irGSEA::irGSEA.density.scatterplot(
          object = GSEA,
          method = "AUCell",
          show.geneset = gsub("_", "-", input$select_geneset),
          reduction = "tsne"
          )
          }
        )
        # VlnPlot
        output$Enri_VlnPlot <- renderPlot(
          {
          irGSEA::irGSEA.halfvlnplot(
          object = GSEA,
          method = "AUCell",
          show.geneset = gsub("_", "-", input$select_geneset)
          )
          }
        )
        }
      )
      }
    )
    }
  )
  }
  ## run shiny
  shinyApp(ui, server)
}
