tabItem(
tabName = "home",
fluidRow(
  column(
    width = 12,
    h2("Welcome to iWGCNA"),
    hr()
  )
),
fluidRow(
  column(
    width = 6,
    uiOutput("background"),
    uiOutput("user")
  ),
  column(
    width = 6,
    uiOutput("provide"),
    uiOutput("highlight"),
    uiOutput("access")
  )
),
fluidRow(
  column(
    width = 12,
    hr(),
    img(
      src = "outline.jpg", 
      style = "display: block; margin: 0 auto; box-shadow: 0 8px 20px rgba(0,0,0,0.4); border-radius: 8px;", 
      width = "100%"
      ),
    )
  )
)