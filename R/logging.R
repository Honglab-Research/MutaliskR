# Logging function.
#
# Author: Andy Jinseok Lee


#' @title Prints log
#'
#' @description This function prints a log message.
#'
#' @param message String value message to print along with log type and date.
#' @param type String value that represents type of this message. 'INFO' by default.
#'
#' @export
PrintLog <- function(message, type = "INFO") {
  print(paste0("[", type, "][", Sys.time(), "] ", message))
}
