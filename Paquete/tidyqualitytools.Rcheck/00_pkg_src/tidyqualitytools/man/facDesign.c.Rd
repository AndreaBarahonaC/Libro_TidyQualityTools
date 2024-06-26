% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0.2 utilsfacDesign.R
\name{facDesign.c}
\alias{facDesign.c}
\title{facDesign Class}
\description{
facDesign Class
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{name}}{name of de facDesign}

\item{\code{factors}}{description}

\item{\code{cube}}{description}

\item{\code{star}}{description}

\item{\code{centerCube}}{description}

\item{\code{centerStar}}{description}

\item{\code{generator}}{description}

\item{\code{response}}{description}

\item{\code{block}}{description}

\item{\code{blockGen}}{description}

\item{\code{runOrder}}{description}

\item{\code{standardOrder}}{description}

\item{\code{desireVal}}{description}

\item{\code{desirability}}{description}

\item{\code{fits}}{description}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-facDesign-nrow}{\code{facDesign.c$nrow()}}
\item \href{#method-facDesign-ncol}{\code{facDesign.c$ncol()}}
\item \href{#method-facDesign-print}{\code{facDesign.c$print()}}
\item \href{#method-facDesign-.clear}{\code{facDesign.c$.clear()}}
\item \href{#method-facDesign-names}{\code{facDesign.c$names()}}
\item \href{#method-facDesign-as.data.frame}{\code{facDesign.c$as.data.frame()}}
\item \href{#method-facDesign-get}{\code{facDesign.c$get()}}
\item \href{#method-facDesign-lows}{\code{facDesign.c$lows()}}
\item \href{#method-facDesign-highs}{\code{facDesign.c$highs()}}
\item \href{#method-facDesign-.nfp}{\code{facDesign.c$.nfp()}}
\item \href{#method-facDesign-identity}{\code{facDesign.c$identity()}}
\item \href{#method-facDesign-summary}{\code{facDesign.c$summary()}}
\item \href{#method-facDesign-.response}{\code{facDesign.c$.response()}}
\item \href{#method-facDesign-effectPlot}{\code{facDesign.c$effectPlot()}}
\item \href{#method-facDesign-lm}{\code{facDesign.c$lm()}}
\item \href{#method-facDesign-desires}{\code{facDesign.c$desires()}}
\item \href{#method-facDesign-set.fits}{\code{facDesign.c$set.fits()}}
\item \href{#method-facDesign-types}{\code{facDesign.c$types()}}
\item \href{#method-facDesign-unit}{\code{facDesign.c$unit()}}
\item \href{#method-facDesign-.star}{\code{facDesign.c$.star()}}
\item \href{#method-facDesign-.blockGen}{\code{facDesign.c$.blockGen()}}
\item \href{#method-facDesign-.block}{\code{facDesign.c$.block()}}
\item \href{#method-facDesign-.centerCube}{\code{facDesign.c$.centerCube()}}
\item \href{#method-facDesign-.centerStar}{\code{facDesign.c$.centerStar()}}
\item \href{#method-facDesign-.generator}{\code{facDesign.c$.generator()}}
\item \href{#method-facDesign-clone}{\code{facDesign.c$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-nrow"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-nrow}{}}}
\subsection{Method \code{nrow()}}{
Get the number of row Design
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$nrow()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-ncol"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-ncol}{}}}
\subsection{Method \code{ncol()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$ncol()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-print"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.clear"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.clear}{}}}
\subsection{Method \code{.clear()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.clear()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-names"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-names}{}}}
\subsection{Method \code{names()}}{
Set the factors names
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$names(value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{value}}{factors names}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-as.data.frame"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-as.data.frame}{}}}
\subsection{Method \code{as.data.frame()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$as.data.frame(row.names = NULL, optional = FALSE, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-get"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-get}{}}}
\subsection{Method \code{get()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$get(i, j)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-lows"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-lows}{}}}
\subsection{Method \code{lows()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$lows(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-highs"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-highs}{}}}
\subsection{Method \code{highs()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$highs(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.nfp"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.nfp}{}}}
\subsection{Method \code{.nfp()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.nfp()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-identity"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-identity}{}}}
\subsection{Method \code{identity()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$identity()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-summary"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-summary}{}}}
\subsection{Method \code{summary()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$summary()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.response"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.response}{}}}
\subsection{Method \code{.response()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.response(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-effectPlot"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-effectPlot}{}}}
\subsection{Method \code{effectPlot()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$effectPlot(
  factors,
  fun = mean,
  response = NULL,
  single = FALSE,
  points = FALSE,
  classic = FALSE,
  axes = TRUE,
  lty,
  xlab,
  ylab,
  main,
  ylim,
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-lm"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-lm}{}}}
\subsection{Method \code{lm()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$lm(formula)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-desires"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-desires}{}}}
\subsection{Method \code{desires()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$desires(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-set.fits"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-set.fits}{}}}
\subsection{Method \code{set.fits()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$set.fits(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-types"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-types}{}}}
\subsection{Method \code{types()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$types(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-unit"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-unit}{}}}
\subsection{Method \code{unit()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$unit(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.star"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.star}{}}}
\subsection{Method \code{.star()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.star(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.blockGen"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.blockGen}{}}}
\subsection{Method \code{.blockGen()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.blockGen(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.block"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.block}{}}}
\subsection{Method \code{.block()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.block(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.centerCube"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.centerCube}{}}}
\subsection{Method \code{.centerCube()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.centerCube(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.centerStar"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.centerStar}{}}}
\subsection{Method \code{.centerStar()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.centerStar(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-.generator"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-.generator}{}}}
\subsection{Method \code{.generator()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$.generator(value)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-facDesign-clone"></a>}}
\if{latex}{\out{\hypertarget{method-facDesign-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{facDesign.c$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
