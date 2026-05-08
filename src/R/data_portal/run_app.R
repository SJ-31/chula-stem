library(shiny)

runApp(
  appDir = getwd(),
  port = 4214,
  launch.browser = TRUE,
  workerId = "",
  quiet = FALSE,
  display.mode = c("auto", "normal", "showcase"),
  test.mode = getOption("shiny.testmode", FALSE)
)
