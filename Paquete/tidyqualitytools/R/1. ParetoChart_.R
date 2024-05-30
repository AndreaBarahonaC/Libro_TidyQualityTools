

############# Nueva Función Pareto (se quita argumento `las`)################
ParetoChart_ <- function (x, weight, showTable = TRUE, showPlot = TRUE, main, col, border, xlab, ylab = "Frequency", percentVec, ...)
{
  varName = deparse(substitute(x))[1]
  corp.col = "#C4B9FF"
  corp.border = "#9E0138"
  if (!is.vector(x) & !is.data.frame(x) & !is.table(x))
    stop("x should be a vector, dataframe or a table")
  if (is.table(x)) {
    xtable = x
  }
  if (is.vector(x)) {
    if (!is.null(names(x)))
      xtable = as.table(x)
    else xtable = table(x)
  }
  if (!missing(weight)) {
    if (!is.numeric(weight))
      stop("weight must be numeric!")
    if (is.null(names(weight)))
      stop("weight is missing names for matching!")
    else {
      if (FALSE %in% (sort(names(weight)) == sort(names(xtable))))
        stop("names of weight and table do not match!")
      else {
        for (i in 1:length(xtable)) {
          xtable[i] = weight[names(weight) == names(xtable)[i]] *
            xtable[i]
        }
      }
    }
  }
  else {
    weight = FALSE
  }
  if (missing(showTable))
    showTable = TRUE
  if (missing(xlab))
    xlab = ""
  if (missing(main))
    main = paste("Pareto Chart for", varName)
  if (missing(col))
    col = corp.col
  if (missing(border))
    border = corp.border
  if (missing(percentVec))
    percentVec = seq(0, 1, by = 0.25)
  call <- match.call(expand.dots = TRUE)
  # Plot
  if (length(xtable) > 1) {
    ylim = c(min(xtable), max(xtable) * 1.025)
    xtable = c(sort(xtable, decreasing = TRUE, na.last = TRUE))
    cumFreq = cumsum(xtable)
    sumFreq = sum(xtable)
    percentage = xtable/sum(xtable) * 100
    cumPerc = cumFreq/sumFreq * 100


    data <- data.frame(Frequency = xtable,
                       Cum.Frequency = cumFreq,
                       Percentage = round(percentage, digits = 2),
                       Cum.Percentage = round(cumPerc, digits = 2))
    tabla <- t(data)

    p <- ggplot(data, aes(x = reorder(names(xtable), -xtable), y = Frequency)) +
      geom_col(aes(fill = "Frequency"), width = 0.7) +
      geom_point(aes(y = Cum.Frequency, color = "Cumulative Percentage"), size = 3) +
      geom_line(aes(y = Cum.Frequency, group = 1, color = "Cumulative Percentage"), line = 0.5) +
      scale_y_continuous(name = ylab,
                         sec.axis = sec_axis(~ . / sum(xtable),
                                             name = "Cumulative Percentage",
                                             labels = percentVec)) +
      scale_x_discrete(name = xlab) +
      scale_color_manual(values = c(border, border)) +
      scale_fill_manual(values = col) +
      theme(legend.position = "none") +
      labs(title = main)+theme(plot.title = element_text(hjust = 0.5,face = "bold"))
  }
  else {
    warning("data should have at least two categories!")
  }
  if(showPlot == TRUE){
    if(showTable == TRUE){
      show(p/tableGrob(tabla))
    }
    else {
      show(p)
    }
  }
  else{
    show(tabla)
  }
}
