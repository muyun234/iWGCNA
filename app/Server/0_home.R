# server

observeEvent(input$home_start, {
  updateTabItems(session, "steps", "upload")
})

observeEvent(input$home_demo, {
  updateTabItems(session, "steps", "upload")
  updateAwesomeCheckbox(session, "demo_on", value = TRUE)
})