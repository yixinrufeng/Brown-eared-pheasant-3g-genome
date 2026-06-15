# ============================================================
# ParFileViewer_R441_safe.R
# 适用于 R 4.4.1 的 fastsimcoal .par 文件绘图脚本
# 主要修复：
# 1) 对空格/tab/多余空格/行尾注释更稳健
# 2) 对 keep、nomig、instbot、absoluteResize 更稳健
# 3) 避免 historical events 中 NA 导致 if(...) 报错
# 4) 如果 .par 某行格式有问题，会直接报告具体行号
# ============================================================

# ------------------------- 0. 安装和载入包 -------------------------

pkgs <- c("plotrix", "diagram", "shape")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(plotrix)
  library(diagram)
  library(shape)
})

# ------------------------- 1. 基本参数 -------------------------

# 不要改
mutRateRescaling <- 1
genTime <- 1
rescalingFactor <- mutRateRescaling * genTime

# 绘图参数，可按需修改
migrMatCol <- "coral"           # gene flow arrows
admixCol <- "blue"              # admixture events
popFusionColor <- "black"       # population divergence/fusion lines
popCol <- "gray95"              # population size rectangles
popBorderCol <- "black"
timeCol <- "tan4"
ageCol <- "dodgerblue2"
textInstbotColor <- "red"

t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  rgb(
    rgb.val[1], rgb.val[2], rgb.val[3],
    max = 255,
    alpha = (100 - percent) * 255 / 100,
    names = name
  )
}

instbotColor <- t_col(textInstbotColor)

propLastsegment <- 0.05
migMatNameProp <- 0.8
migMatLineLength <- 0.3
timeProp <- 0.6
maxRadius <- 1 / 40
minRadius <- maxRadius / 3
arrowLength <- 0.2
timeOffset <- 0.25
migrOffset0 <- 0.05
curvedArrowLTY <- 1
drawLogPopSize <- FALSE
plotMigrRates <- TRUE
plotNmValues <- FALSE
migrRateTextSizeCEX <- 0.5
instBotAbsSize <- TRUE
maxNumSpecialEvents <- 50

pdf.x.size <- 10
pdf.y.size <- 8

# ------------------------- 2. 文件设置 -------------------------

args <- commandArgs(TRUE)

if (length(args) >= 1 && !is.na(args[1]) && nzchar(args[1])) {

  # 用命令行运行：
  # Rscript ParFileViewer_R441_safe.R your_model.par
  parFileName <- args[1]
  printPDF <- TRUE
  usePopLabels <- FALSE
  fixMaxTime <- FALSE
  externalCall <- TRUE

} else {

  # 在 RStudio 里运行时，主要改这里
  printPDF <- FALSE

  # 改成你的 .par 文件所在文件夹
  setwd("D:/temp/新建文件夹 (2)")

  # 改成你的 .par 文件名
  parFileName <- "contr.par"

  # 如果有两个种群，可以这样写；如果数量不对，脚本会自动改用 0,1,2...
  popLabels <- c("Pop1", "Pop2")
  usePopLabels <- TRUE

  # 是否限制最大绘图时间
  fixMaxTime <- TRUE
  maxTimeToPlot <- 8000

  externalCall <- FALSE
}

if (!exists("popLabels")) {
  usePopLabels <- FALSE
}

# ------------------------- 3. 稳健读取函数 -------------------------

clean_line <- function(x) {
  x <- gsub("#.*$", "", x)
  x <- gsub("//.*$", "", x)
  x <- trimws(x)
  x <- gsub("\t", " ", x)
  x <- gsub(" +", " ", x)
  x
}

tokens_of <- function(x) {
  x <- clean_line(x)
  if (is.na(x) || !nzchar(x)) return(character(0))
  unlist(strsplit(x, " ", fixed = TRUE))
}

replace_keep_tokens <- function(x) {
  x[grepl("^keep$", x, ignore.case = TRUE)] <- "-9999"
  x
}

stop_line <- function(message, line_no, raw_line) {
  stop(
    paste0(
      "\n", message, "\n",
      "par 文件行号: ", line_no, "\n",
      "原始内容: ", raw_line, "\n"
    ),
    call. = FALSE
  )
}

first_number <- function(parFile, line_no, label = "number") {
  if (line_no < 1 || line_no > length(parFile)) {
    stop(paste0("读取 ", label, " 时超出 par 文件行数。"), call. = FALSE)
  }

  tok <- replace_keep_tokens(tokens_of(parFile[line_no]))

  if (length(tok) < 1) {
    stop_line(paste0("读取 ", label, " 失败：该行为空。"), line_no, parFile[line_no])
  }

  val <- suppressWarnings(as.numeric(tok[1]))

  if (is.na(val)) {
    stop_line(
      paste0("读取 ", label, " 失败：第一列不是数字。"),
      line_no,
      parFile[line_no]
    )
  }

  val
}

numeric_fields <- function(parFile, line_no, min_n = 1, label = "numeric fields") {
  if (line_no < 1 || line_no > length(parFile)) {
    stop(paste0("读取 ", label, " 时超出 par 文件行数。"), call. = FALSE)
  }

  raw <- parFile[line_no]
  tok <- replace_keep_tokens(tokens_of(raw))

  if (length(tok) < min_n) {
    stop_line(
      paste0("读取 ", label, " 失败：字段数量不足，至少需要 ", min_n, " 列。"),
      line_no,
      raw
    )
  }

  nums <- suppressWarnings(as.numeric(tok[seq_len(min_n)]))

  if (any(is.na(nums))) {
    stop_line(
      paste0(
        "读取 ", label, " 失败：前 ", min_n,
        " 列中有字段不能转换为数字。\n",
        "解析字段: ", paste(tok, collapse = " | ")
      ),
      line_no,
      raw
    )
  }

  nums
}

getNumbers <- function(start, parFile, numSamples, label = "numbers") {
  out <- numeric(numSamples)

  for (i in seq_len(numSamples)) {
    line_no <- start + i
    out[i] <- first_number(
      parFile = parFile,
      line_no = line_no,
      label = paste0(label, " line ", i)
    )
  }

  out
}

readSampleSizesTimesAndInbreedingLevel <- function(start, parFile, numSamples) {
  sampSize <- numeric(numSamples)
  sampTime <- numeric(numSamples)
  inbreeding <- numeric(numSamples)

  for (i in seq_len(numSamples)) {
    line_no <- start + i
    raw <- parFile[line_no]
    tok <- replace_keep_tokens(tokens_of(raw))

    if (length(tok) < 1) {
      stop_line("读取 sample size 失败：该行为空。", line_no, raw)
    }

    nums <- suppressWarnings(as.numeric(tok))

    if (is.na(nums[1])) {
      stop_line("读取 sample size 失败：第一列不是数字。", line_no, raw)
    }

    sampSize[i] <- nums[1]
    sampTime[i] <- ifelse(length(nums) >= 2 && !is.na(nums[2]), nums[2], 0)
    inbreeding[i] <- ifelse(length(nums) >= 3 && !is.na(nums[3]), nums[3], 0)
  }

  list(ss = sampSize, st = sampTime, inb = inbreeding)
}

readMigMat <- function(start, parFile, numSamples, mat_no = 1) {
  migs <- matrix(0, nrow = numSamples, ncol = numSamples)

  for (i in seq_len(numSamples)) {
    line_no <- start + i
    row_vals <- numeric_fields(
      parFile = parFile,
      line_no = line_no,
      min_n = numSamples,
      label = paste0("migration matrix ", mat_no, " row ", i)
    )
    migs[i, ] <- row_vals
  }

  rownames(migs) <- seq_len(numSamples)
  colnames(migs) <- seq_len(numSamples)
  migs
}

# ------------------------- 4. 读取 par 文件 -------------------------

if (!file.exists(parFileName)) {
  stop(paste0("找不到 par 文件：", parFileName), call. = FALSE)
}

parFile <- scan(
  parFileName,
  character(0),
  sep = "\n",
  strip.white = FALSE,
  quiet = TRUE
)

if (length(parFile) < 2) {
  stop("par 文件太短，无法读取。", call. = FALSE)
}

if (fixMaxTime) {
  suffixZoom <- paste("_zoom_", maxTimeToPlot, sep = "")
} else {
  suffixZoom <- ""
}

if (genTime == 1) {
  suffixGen <- "_gens"
} else {
  suffixGen <- "_years"
}

if (externalCall) suffixGen <- ""

pdfFileName <- paste(parFileName, suffixZoom, suffixGen, ".pdf", sep = "")

# line 2: number of samples/demes
numSamples <- first_number(parFile, 2, "number of samples")
numSamples <- as.integer(numSamples)
numDemes <- numSamples

if (exists("popLabels") && length(popLabels) >= numSamples) {
  usePopLabels <- TRUE
} else {
  usePopLabels <- FALSE
}

# ------------------------- 5. 读取初始种群大小 -------------------------

numPopSizes <- rep(1, numDemes)
popSizes <- array(0, dim = c(maxNumSpecialEvents, numDemes))
popSizeTimeChange <- array(0, dim = c(maxNumSpecialEvents, numDemes))

isGrowing <- rep(FALSE, numDemes)
curGrowthRates <- rep(0, numDemes)
lastGrowthStartTime <- rep(-1, numDemes)

numGrowths <- rep(0, numDemes)
startGrowthTimes <- array(0, dim = c(maxNumSpecialEvents, numDemes))
endGrowthTimes <- array(0, dim = c(maxNumSpecialEvents, numDemes))
startGrowthSize <- array(0, dim = c(maxNumSpecialEvents, numDemes))
endGrowthSize <- array(0, dim = c(maxNumSpecialEvents, numDemes))

start <- 3

popSizes[1, ] <- getNumbers(
  start = start,
  parFile = parFile,
  numSamples = numSamples,
  label = "initial population size"
)

popSizes[1, ] <- round(popSizes[1, ] * rescalingFactor, digits = 0)
iniPopSizes <- popSizes[1, ]

# ------------------------- 6. 读取样本大小、样本时间、近交系数 -------------------------

start <- start + numSamples + 1

sampSizesStats <- readSampleSizesTimesAndInbreedingLevel(
  start = start,
  parFile = parFile,
  numSamples = numSamples
)

sampSizesStats$st <- round(sampSizesStats$st * genTime, digits = 0)

sampSizes <- sampSizesStats$ss
sampTimes <- sampSizesStats$st
inbrCoeff <- sampSizesStats$inb

popSizeTimeChange[1, ] <- sampTimes

# ------------------------- 7. 读取 growth rates -------------------------

start <- start + numSamples + 1

growthRatesInitial <- getNumbers(
  start = start,
  parFile = parFile,
  numSamples = numSamples,
  label = "growth rate"
)

growthRates <- growthRatesInitial
curGrowthRates <- growthRatesInitial

for (i in seq_len(numDemes)) {
  if (!is.na(curGrowthRates[i]) && curGrowthRates[i] != 0) {
    isGrowing[i] <- TRUE
    lastGrowthStartTime[i] <- 0
    numGrowths[i] <- numGrowths[i] + 1
    startGrowthTimes[numGrowths[i], i] <- sampTimes[i]
    startGrowthSize[numGrowths[i], i] <- iniPopSizes[i]
  }
}

# ------------------------- 8. 读取 migration matrices -------------------------

start <- start + numSamples + 1

numMigMat <- first_number(parFile, start + 1, "number of migration matrices")
numMigMat <- as.integer(numMigMat)

start <- start + 2

migMats <- list()

if (!is.na(numMigMat) && numMigMat > 0) {
  for (i in seq_len(numMigMat)) {
    migMats[[i]] <- readMigMat(
      start = start,
      parFile = parFile,
      numSamples = numSamples,
      mat_no = i
    )
    start <- start + numSamples + 1
  }
}

# ------------------------- 9. 读取 historical events -------------------------

start <- start + 1

numHistEvents <- first_number(parFile, start, "number of historical events")
numHistEvents <- as.integer(numHistEvents)

timesPos <- 0
migrMatPos <- numSamples + 1

last.he.time <- 0
histEvents <- NULL

if (!is.na(numHistEvents) && numHistEvents > 0) {

  histEvents <- matrix(NA_real_, nrow = numHistEvents, ncol = 9)

  for (i in seq_len(numHistEvents)) {

    start <- start + 1
    raw_line <- parFile[start]
    fields <- replace_keep_tokens(tokens_of(raw_line))

    if (length(fields) < 7) {
      stop_line(
        paste0(
          "Historical event 读取失败：字段数量不足。\n",
          "事件编号: ", i, "\n",
          "至少需要 7 列：time source sink migrants resize newGrowth newMigMat"
        ),
        start,
        raw_line
      )
    }

    nums7 <- suppressWarnings(as.numeric(fields[1:7]))

    if (any(is.na(nums7))) {
      stop_line(
        paste0(
          "Historical event 读取失败：前 7 列必须是数字或 keep。\n",
          "事件编号: ", i, "\n",
          "解析字段: ", paste(fields, collapse = " | "), "\n",
          "标准格式通常是：\n",
          "time source sink migrants resize newGrowth newMigMat\n",
          "如果使用关键词 nomig、instbot、absoluteResize，建议放在前 7 列之后。"
        ),
        start,
        raw_line
      )
    }

    curHE <- numeric(9)
    curHE[1:7] <- nums7

    nomig <- any(grepl("^nomig$", fields, ignore.case = TRUE))
    instbot <- any(grepl("^instbot$", fields, ignore.case = TRUE))
    absResize <- any(grepl("^absoluteresize$", fields, ignore.case = TRUE))

    if (nomig) curHE[7] <- -1
    curHE[8] <- ifelse(instbot, 1, 0)
    curHE[9] <- ifelse(absResize, 1, 0)

    curHE[1] <- round(curHE[1] * genTime, digits = 0)

    histEvents[i, ] <- curHE

    if (!is.na(curHE[1]) && !is.na(last.he.time) && curHE[1] > last.he.time) {
      last.he.time <- curHE[1]
    }
  }

  colnames(histEvents) <- c(
    "time", "source", "sink", "migr", "resize",
    "growth", "migmat", "instBot", "absResize"
  )

  rownames(histEvents) <- seq_len(numHistEvents)

  # 按时间排序，避免 par 文件中 historical events 顺序不规范
  histEvents <- histEvents[order(histEvents[, "time"], decreasing = FALSE), , drop = FALSE]
}

endReadParFile <- start

# ------------------------- 10. 绘图辅助函数 -------------------------

fullHeadArrow <- function(x0, y0, x1, y1, length, angle, color = "black", weight = 1) {
  arrows(x0, y0, x1, y1, length, angle, code = 2, lty = 1, col = color, lwd = weight)
  arrows(x0, y0, x1, y1, length, angle * 0.80, code = 2, lty = 1, col = color, lwd = weight)
  arrows(x0, y0, x1, y1, length, angle * 0.60, code = 2, lty = 1, col = color, lwd = weight)
  arrows(x0, y0, x1, y1, length, angle * 0.40, code = 2, lty = 1, col = color, lwd = weight)
  arrows(x0, y0, x1, y1, length, angle * 0.20, code = 2, lty = 1, col = color, lwd = weight)
  arrows(x0, y0, x1, y1, length, angle * 0.10, code = 2, lty = 1, col = color, lwd = weight)
}

drawTriangle <- function(growth, x, y, size, aspRatio, color) {
  if (growth > 0) {
    x0 <- x
    y0 <- y
    x1 <- x - size / 2
    y1 <- y + size / 2 * aspRatio
    x2 <- x + size / 2
    y2 <- y + size / 2 * aspRatio
  } else {
    x0 <- x - size / 2
    y0 <- y
    x1 <- x
    y1 <- y + size / 2 * aspRatio
    x2 <- x + size / 2
    y2 <- y
  }

  polygon(c(x0, x1, x2, x0), c(y0, y1, y2, y0), col = color)
  y + size / 2 * aspRatio
}

interpolRadius <- function(curSize, minSize, maxSize, minRadius, maxRadius, logScale) {
  if (!is.finite(curSize) || !is.finite(minSize) || !is.finite(maxSize)) {
    return((minRadius + maxRadius) / 2)
  }

  if (logScale) {
    minSize <- ifelse(minSize > 0, log10(minSize), 0)
    maxSize <- ifelse(maxSize > 0, log10(maxSize), 0)
    curSize <- ifelse(curSize > 0, log10(curSize), 0)
  }

  if (maxSize == minSize) {
    return((minRadius + maxRadius) / 2)
  }

  minRadius + (curSize - minSize) * (maxRadius - minRadius) / (maxSize - minSize)
}

getPopSize <- function(time, sink, allPopSizes, timePopSizes, numChanges) {
  if (sink < 1 || sink > ncol(allPopSizes)) return(0)

  if (numChanges[sink] <= 1) {
    return(allPopSizes[1, sink])
  }

  for (i in 2:numChanges[sink]) {
    if (timePopSizes[i, sink] > time) {
      return(allPopSizes[i - 1, sink])
    }
  }

  allPopSizes[numChanges[sink], sink]
}

# ------------------------- 11. 计算种群大小变化历史 -------------------------

maxPopSize <- max(iniPopSizes, na.rm = TRUE)
minPopSize <- min(iniPopSizes, na.rm = TRUE)

numInstBot <- rep(0, numDemes)
sizeInstBot <- array(NA_real_, dim = c(maxNumSpecialEvents, numDemes))
timeInstBot <- array(0, dim = c(maxNumSpecialEvents, numDemes))

numMigEvents <- rep(0, numDemes)
sizeMigr <- array(0, dim = c(maxNumSpecialEvents, numDemes))
timeMigEvent <- array(0, dim = c(maxNumSpecialEvents, numDemes))
sourceMigEvent <- array(0, dim = c(maxNumSpecialEvents, numDemes))

curMigMatNum <- 0
numMigmatChanges <- 0
migMatNumbers <- rep(0, maxNumSpecialEvents)
timeMigMatChanges <- rep(0, maxNumSpecialEvents)

activeAtEnd <- rep(TRUE, numDemes)
curTime <- 0

if (!is.null(histEvents) && nrow(histEvents) > 0) {

  for (i in seq_len(nrow(histEvents))) {

    he <- histEvents[i, ]

    curTime <- he["time"]
    source <- he["source"] + 1
    sink <- he["sink"] + 1
    migr <- he["migr"]
    resize <- he["resize"]
    growth <- he["growth"]
    migrMat <- he["migmat"]
    isInstBot <- he["instBot"]
    isAbsResize <- he["absResize"]

    if (source < 1 || source > numDemes || sink < 1 || sink > numDemes) {
      stop(
        paste0(
          "\nHistorical event 中 source/sink 超出种群编号范围。\n",
          "事件编号: ", i, "\n",
          "source = ", source - 1, ", sink = ", sink - 1, "\n",
          "注意：fastsimcoal 的 deme 编号通常从 0 开始，最大为 numSamples - 1。\n"
        ),
        call. = FALSE
      )
    }

    if (growth == -9999) {
      growth <- curGrowthRates[sink]
    }

    if (migrMat == -9999) {
      migrMat <- curMigMatNum
    }

    # 记录 migration matrix changes
    if (!is.na(migrMat) && migrMat != curMigMatNum) {
      numMigmatChanges <- numMigmatChanges + 1
      if (numMigmatChanges > maxNumSpecialEvents) {
        stop("migration matrix changes 太多，请增大 maxNumSpecialEvents。", call. = FALSE)
      }
      migMatNumbers[numMigmatChanges] <- migrMat
      timeMigMatChanges[numMigmatChanges] <- curTime
      curMigMatNum <- migrMat
    }

    # 处理 growth
    if (growth != -9999) {

      if (isGrowing[sink] && growth == 0) {

        lastPopSize <- popSizes[numPopSizes[sink], sink]
        deltaT <- curTime - popSizeTimeChange[numPopSizes[sink], sink]
        newpopSize <- round(lastPopSize * exp(curGrowthRates[sink] * deltaT / genTime))

        curGrowthRates[sink] <- 0
        isGrowing[sink] <- FALSE

        numPopSizes[sink] <- numPopSizes[sink] + 1
        popSizeTimeChange[numPopSizes[sink], sink] <- curTime
        popSizes[numPopSizes[sink], sink] <- newpopSize

      } else if (!isGrowing[sink] && growth != 0) {

        isGrowing[sink] <- TRUE
        curGrowthRates[sink] <- growth

        numPopSizes[sink] <- numPopSizes[sink] + 1
        popSizeTimeChange[numPopSizes[sink], sink] <- curTime
        popSizes[numPopSizes[sink], sink] <- popSizes[numPopSizes[sink] - 1, sink]

      } else if (isGrowing[sink] && growth != 0) {

        lastPopSize <- popSizes[numPopSizes[sink], sink]
        deltaT <- curTime - popSizeTimeChange[numPopSizes[sink], sink]
        newpopSize <- round(lastPopSize * exp(curGrowthRates[sink] * deltaT / genTime))

        numPopSizes[sink] <- numPopSizes[sink] + 1
        popSizeTimeChange[numPopSizes[sink], sink] <- curTime
        popSizes[numPopSizes[sink], sink] <- newpopSize

        curGrowthRates[sink] <- growth
      }
    }

    # 记录 size change
    if (!is.na(resize) && resize != 1 && isInstBot == 0) {

      oldPopSize <- popSizes[numPopSizes[sink], sink]

      numPopSizes[sink] <- numPopSizes[sink] + 1
      popSizeTimeChange[numPopSizes[sink], sink] <- curTime
      popSizes[numPopSizes[sink], sink] <- oldPopSize

      numPopSizes[sink] <- numPopSizes[sink] + 1
      popSizeTimeChange[numPopSizes[sink], sink] <- curTime

      if (isAbsResize == 0) {
        popSizes[numPopSizes[sink], sink] <- round(oldPopSize * resize)
      } else {
        popSizes[numPopSizes[sink], sink] <- round(resize)
      }
    }

    # 处理 population fusion
    if (!is.na(migr) && migr >= 1 && source != sink) {

      if (!isGrowing[source]) {
        numPopSizes[source] <- numPopSizes[source] + 1
        popSizeTimeChange[numPopSizes[source], source] <- curTime
        popSizes[numPopSizes[source], source] <- popSizes[numPopSizes[source] - 1, source]
      } else {
        endGrowthTimes[numGrowths[source], source] <- curTime
        lastPopSize <- popSizes[numPopSizes[source], source]
        deltaT <- curTime - startGrowthTimes[numGrowths[source], source]
        newpopSize <- round(lastPopSize * exp(curGrowthRates[source] * deltaT / genTime))
        endGrowthSize[numGrowths[source], source] <- newpopSize
      }

      activeAtEnd[source] <- FALSE
    }

    # instantaneous bottlenecks
    if (isInstBot == 1) {
      numInstBot[sink] <- numInstBot[sink] + 1
      if (numInstBot[sink] > maxNumSpecialEvents) {
        stop("instantaneous bottleneck 事件太多，请增大 maxNumSpecialEvents。", call. = FALSE)
      }
      timeInstBot[numInstBot[sink], sink] <- curTime
      sizeInstBot[numInstBot[sink], sink] <- 1 / resize
    }

    # admixture / pulse migration events
    if (!is.na(migr) && migr > 0 && migr < 1) {
      numMigEvents[sink] <- numMigEvents[sink] + 1
      if (numMigEvents[sink] > maxNumSpecialEvents) {
        stop("migration events 太多，请增大 maxNumSpecialEvents。", call. = FALSE)
      }
      timeMigEvent[numMigEvents[sink], sink] <- curTime
      sizeMigr[numMigEvents[sink], sink] <- migr
      sourceMigEvent[numMigEvents[sink], sink] <- source
    }

    curSize <- popSizes[numPopSizes[sink], sink]

    if (is.finite(curSize)) {
      maxPopSize <- max(maxPopSize, curSize, na.rm = TRUE)
      minPopSize <- min(minPopSize, curSize, na.rm = TRUE)
    }

    if (curSize == 0) {
      message("Deme ", sink, " reaches size zero at time ", curTime)
    }

    if (is.infinite(curSize)) {
      message("Deme ", sink, " reaches infinite size at time ", curTime)
    }
  }

  # 将最终仍然存在的 deme 延伸到图形上方
  top_time <- ifelse(curTime > 0, curTime + 0.1 * curTime, 1)
  active_demes <- which(activeAtEnd)

  if (length(active_demes) == 0) active_demes <- sink

  for (d in active_demes) {
    numPopSizes[d] <- numPopSizes[d] + 1
    popSizeTimeChange[numPopSizes[d], d] <- top_time
    popSizes[numPopSizes[d], d] <- popSizes[numPopSizes[d] - 1, d]
  }
}

last.he.time <- as.numeric(last.he.time)

maxPopSizeTimeChange <- suppressWarnings(
  max(c(popSizeTimeChange, timeInstBot), na.rm = TRUE)
)

if (!is.finite(maxPopSizeTimeChange)) maxPopSizeTimeChange <- 0

if (last.he.time != 0) {
  yTimeLimit <- max(last.he.time, maxPopSizeTimeChange) * (1 + propLastsegment)
} else {
  yTimeLimit <- max(maxPopSizeTimeChange, max(sampTimes, na.rm = TRUE), 1)
}

if (fixMaxTime) {
  if (yTimeLimit > maxTimeToPlot) yTimeLimit <- maxTimeToPlot
} else {
  maxTimeToPlot <- yTimeLimit
}

if (!is.finite(yTimeLimit) || yTimeLimit <= 0) yTimeLimit <- 1

# ------------------------- 12. 开始绘图 -------------------------

if (printPDF) {
  pdf(pdfFileName, width = pdf.x.size, height = pdf.y.size)
}

par(xpd = FALSE, mar = c(8, 6, 3, 2))

maxRadius <- maxRadius * (numSamples + 2)
minRadius <- minRadius * (numSamples + 2)

title <- parFileName

ylabel <- ifelse(genTime == 1, "time (gen)", "time (years)")

plot(
  x = seq_len(numSamples),
  type = "n",
  xlab = "",
  ylab = "",
  xlim = c(-0.5, numSamples + 1.5),
  ylim = c(0, yTimeLimit),
  main = title,
  xaxt = "n",
  cex.main = 0.8,
  cex.axis = 0.9,
  cex.lab = 0.9,
  las = 2
)

mtext(side = 2, ylabel, line = 4)

if (usePopLabels) {
  axis(
    side = 1,
    labels = c("Times", popLabels[seq_len(numSamples)], " \nMig Mat"),
    at = 0:(numSamples + 1),
    cex.axis = 0.8
  )
} else {
  axis(
    side = 1,
    labels = c("Times", 0:(numSamples - 1), " \nMig Mat"),
    at = 0:(numSamples + 1),
    cex.axis = 0.8
  )
}

mtext(paste("n =", sampSizes), side = 1, line = 1.8, at = seq_len(numDemes), cex = 0.7)
mtext(paste("2N =", iniPopSizes), side = 1, line = 2.5, at = seq_len(numDemes), cex = 0.7)
mtext(
  paste("FIS =", format(inbrCoeff, scientific = FALSE, digits = 2)),
  side = 1,
  line = 3.2,
  at = seq_len(numDemes),
  cex = 0.7
)

w <- par("pin")[1] / diff(par("usr")[1:2])
h <- par("pin")[2] / diff(par("usr")[3:4])
aspRatio <- w / h

slide <- timeOffset

# ------------------------- 13. 画 population size changes -------------------------

for (i in seq_len(numDemes)) {

  if (numPopSizes[i] >= 2) {

    for (j in 2:numPopSizes[i]) {

      if (popSizes[j, i] > 0) {

        curIniRadius <- interpolRadius(
          popSizes[j - 1, i],
          minPopSize,
          maxPopSize,
          minRadius,
          maxRadius,
          drawLogPopSize
        )

        curEndRadius <- interpolRadius(
          popSizes[j, i],
          minPopSize,
          maxPopSize,
          minRadius,
          maxRadius,
          drawLogPopSize
        )

        bottomLeftX <- i - curIniRadius
        bottomRightX <- i + curIniRadius
        topLeftX <- i - curEndRadius
        topRightX <- i + curEndRadius

        bottomLeftY <- popSizeTimeChange[j - 1, i]
        bottomRightY <- popSizeTimeChange[j - 1, i]
        topLeftY <- popSizeTimeChange[j, i]
        topRightY <- popSizeTimeChange[j, i]

        polygon(
          c(bottomLeftX, topLeftX, topRightX, bottomRightX),
          c(bottomLeftY, topLeftY, topRightY, bottomRightY),
          col = popCol,
          border = popBorderCol
        )

        text(
          timesPos + slide,
          bottomLeftY,
          labels = bottomLeftY,
          cex = timeProp,
          col = timeCol
        )

        slide <- -slide
      }
    }
  }
}

# ------------------------- 14. 画第一套 migration matrix -------------------------

curMigMatNum <- 0
curvature <- 0.0075 * max(last.he.time, yTimeLimit)
migrOffset <- migrOffset0

if (numMigMat > 0 && length(migMats) > 0) {

  if (numMigMat == 1 || numMigmatChanges == 0) {
    time2DrawArrows <- yTimeLimit / 2
  } else {
    time2DrawArrows <- timeMigMatChanges[1] / 2
  }

  if (time2DrawArrows != 0) {

    text(
      migrMatPos - migrOffset,
      time2DrawArrows,
      labels = 0,
      cex = migMatNameProp,
      col = migrMatCol
    )

    curMigMat <- migMats[[1]]

    for (sink in seq_len(numSamples)) {
      for (sourc in seq_len(numSamples)) {

        if (sink != sourc && curMigMat[sourc, sink] > 0) {

          differ <- sourc - sink

          curvedarrow(
            from = c(sourc, time2DrawArrows),
            to = c(sink, time2DrawArrows),
            curve = -curvature * (abs(differ) * 0.55 ^ abs(differ)),
            arr.adj = 1,
            arr.pos = 0.5,
            arr.type = "triangle",
            arr.col = migrMatCol,
            lwd = 1,
            lty = curvedArrowLTY,
            lcol = migrMatCol,
            arr.length = arrowLength
          )

          if (plotMigrRates) {

            if (plotNmValues) {
              valueToPlot <- round(curMigMat[sourc, sink] * iniPopSizes[sourc], digits = 2)
            } else {
              valueToPlot <- format(curMigMat[sourc, sink], digits = 2)
            }

            if (sourc > sink) {
              xPosText <- sourc - 0.2 * abs(differ)
              yPosText <- time2DrawArrows + aspRatio * 0.1
            } else {
              xPosText <- sourc + 0.2 * abs(differ)
              yPosText <- time2DrawArrows - aspRatio * 0.1
            }

            text(
              xPosText,
              yPosText,
              labels = valueToPlot,
              cex = migrRateTextSizeCEX,
              col = migrMatCol
            )
          }
        }
      }
    }
  }
}

# ------------------------- 15. 画 instantaneous bottlenecks -------------------------

if (max(numInstBot) > 0) {

  maxBotSize <- max(sizeInstBot, na.rm = TRUE)
  minBotSize <- min(sizeInstBot, na.rm = TRUE)

  if (maxBotSize == minBotSize) maxBotSize <- 4 * maxBotSize

  for (i in seq_len(numDemes)) {

    if (numInstBot[i] > 0) {

      for (j in seq_len(numInstBot[i])) {

        curRadius <- interpolRadius(
          maxBotSize - sizeInstBot[j, i],
          minBotSize,
          maxBotSize,
          minRadius,
          maxRadius,
          drawLogPopSize
        ) / 2

        curPopSize <- getPopSize(
          timeInstBot[j, i],
          i,
          popSizes,
          popSizeTimeChange,
          numPopSizes
        )

        if (curPopSize > 0) {
          draw.circle(
            i,
            timeInstBot[j, i],
            radius = curRadius,
            col = instbotColor,
            border = instbotColor
          )
        }

        if (instBotAbsSize) {
          botText <- sizeInstBot[j, i]
        } else {
          botText <- 1.0 / sizeInstBot[j, i]
        }

        text(
          i,
          timeInstBot[j, i] * 1.005,
          labels = format(botText, digits = 1),
          cex = timeProp,
          col = textInstbotColor,
          pos = 3
        )

        curPopRadius <- interpolRadius(
          curPopSize,
          minPopSize,
          maxPopSize,
          minRadius,
          maxRadius,
          drawLogPopSize
        )

        text(
          i + curPopRadius,
          timeInstBot[j, i],
          labels = format(timeInstBot[j, i], digits = 1),
          cex = timeProp,
          col = textInstbotColor,
          pos = 4
        )
      }
    }
  }
}

# ------------------------- 16. 画 admixture / pulse migration events -------------------------

for (i in seq_len(numDemes)) {

  if (numMigEvents[i] > 0) {

    for (j in seq_len(numMigEvents[i])) {

      curRadius <- interpolRadius(
        getPopSize(timeMigEvent[j, i], i, popSizes, popSizeTimeChange, numPopSizes),
        minPopSize,
        maxPopSize,
        minRadius,
        maxRadius,
        drawLogPopSize
      )

      if (sourceMigEvent[j, i] > i) {
        target <- i + curRadius
      } else {
        target <- i - curRadius
      }

      segments(
        sourceMigEvent[j, i],
        timeMigEvent[j, i],
        target,
        timeMigEvent[j, i],
        col = admixCol,
        lty = 2
      )

      if (i > sourceMigEvent[j, i]) {
        posArrowX <- target - 0.15
        posTextX <- sourceMigEvent[j, i] + 0.5
      } else {
        posArrowX <- target + 0.15
        posTextX <- sourceMigEvent[j, i] - 0.5
      }

      fullHeadArrow(
        posArrowX,
        timeMigEvent[j, i],
        target,
        timeMigEvent[j, i],
        length = 0.15,
        angle = 20,
        color = admixCol
      )

      text(
        timesPos + slide,
        timeMigEvent[j, i],
        labels = timeMigEvent[j, i],
        cex = timeProp,
        col = admixCol
      )

      slide <- -slide

      text(
        posTextX,
        timeMigEvent[j, i],
        labels = format(sizeMigr[j, i], digits = 3),
        cex = timeProp,
        col = admixCol,
        pos = 3
      )
    }
  }
}

# ------------------------- 17. 画 migration matrix changes -------------------------

if (numMigmatChanges > 0 && length(migMats) > 0) {

  for (i in seq_len(numMigmatChanges)) {

    if (i < numMigmatChanges) {
      time2DrawArrows <- (timeMigMatChanges[i + 1] + timeMigMatChanges[i]) / 2
    } else {
      time2DrawArrows <- (timeMigMatChanges[i] + yTimeLimit) / 2
    }

    curMigMatNumber <- migMatNumbers[i]
    mig_index <- curMigMatNumber + 1

    if (mig_index >= 1 && mig_index <= length(migMats)) {

      curMigMat <- migMats[[mig_index]]
      migrOffset <- -migrOffset

      text(
        migrMatPos - migrOffset,
        time2DrawArrows,
        labels = curMigMatNumber,
        cex = migMatNameProp,
        col = migrMatCol
      )

      for (sink in seq_len(numSamples)) {
        for (sourc in seq_len(numSamples)) {

          if (sink != sourc && curMigMat[sourc, sink] > 0) {

            differ <- sourc - sink

            curvedarrow(
              from = c(sourc, time2DrawArrows),
              to = c(sink, time2DrawArrows),
              curve = curvature * (abs(differ) * 0.55 ^ abs(differ)),
              arr.adj = 1,
              arr.pos = 0.5,
              arr.type = "triangle",
              arr.col = migrMatCol,
              lwd = 1,
              lty = curvedArrowLTY,
              lcol = migrMatCol,
              arr.length = arrowLength
            )

            if (plotMigrRates) {

              if (plotNmValues) {
                curSize <- getPopSize(
                  time2DrawArrows,
                  sourc,
                  popSizes,
                  popSizeTimeChange,
                  numPopSizes
                )
                valueToPlot <- round(curMigMat[sourc, sink] * curSize, digits = 2)
              } else {
                valueToPlot <- format(curMigMat[sourc, sink], digits = 2)
              }

              if (sourc > sink) {
                xPosText <- sourc - 0.2 * abs(differ)
                yPosText <- time2DrawArrows + aspRatio * 0.1
              } else {
                xPosText <- sourc + 0.2 * abs(differ)
                yPosText <- time2DrawArrows - aspRatio * 0.1
              }

              text(
                xPosText,
                yPosText,
                labels = valueToPlot,
                cex = migrRateTextSizeCEX,
                col = migrMatCol
              )
            }
          }
        }
      }

      segments(
        migrMatPos - migMatLineLength / 2,
        timeMigMatChanges[i],
        migrMatPos + migMatLineLength / 2,
        timeMigMatChanges[i],
        lty = 3,
        col = migrMatCol
      )
    } else {
      warning(
        paste0(
          "Historical event 使用了 migration matrix 编号 ",
          curMigMatNumber,
          "，但实际只读取到 ",
          length(migMats),
          " 个 migration matrix。"
        )
      )
    }
  }
}

# ------------------------- 18. 画 population fusion lines -------------------------

lastTime <- 0
activePops <- seq_len(numSamples)
numActivePops <- numSamples
lastSink <- -1
removedPops <- integer(0)

if (!is.null(histEvents) && nrow(histEvents) > 0) {

  for (i in seq_len(nrow(histEvents))) {

    he <- histEvents[i, ]

    he.time <- he["time"]
    he.source <- he["source"] + 1
    he.sink <- he["sink"] + 1
    he.migr <- he["migr"]
    he.resize <- he["resize"]
    he.growth <- he["growth"]
    he.migrMat <- he["migmat"]
    he.instbot <- he["instBot"]

    if (he.growth == -9999) {
      he.growth <- growthRates[he.sink]
    }

    if (he.migrMat == -9999) {
      he.migrMat <- curMigMatNum
    }

    if (!is.na(he.migr) && he.migr >= 1 && he.sink != he.source) {

      removedPops <- unique(c(removedPops, he.source))

      curSize <- getPopSize(
        he.time,
        he.sink,
        popSizes,
        popSizeTimeChange,
        numPopSizes
      )

      curRadius <- interpolRadius(
        curSize,
        minPopSize,
        maxPopSize,
        minRadius,
        maxRadius,
        drawLogPopSize
      )

      numActivePops <- numActivePops - 1
      activePops <- setdiff(seq_len(numSamples), removedPops)

      if (he.source > he.sink) {
        target <- he.sink + curRadius
      } else {
        target <- he.sink - curRadius
      }

      segments(
        he.source,
        he.time,
        target,
        he.time,
        col = popFusionColor,
        lty = 1
      )

      text(
        timesPos + slide,
        he.time,
        labels = he.time,
        cex = timeProp,
        col = popFusionColor
      )

      slide <- -slide
    }

    lastSink <- he.sink
  }
}

# 画最终仍然存在的分支
active_demes <- which(activeAtEnd)
if (length(active_demes) == 0 && lastSink >= 1) active_demes <- lastSink
if (length(active_demes) == 0) active_demes <- seq_len(numDemes)

for (d in active_demes) {
  if (popSizeTimeChange[numPopSizes[d], d] < maxTimeToPlot) {
    segments(
      d,
      popSizeTimeChange[numPopSizes[d], d],
      d,
      yTimeLimit
    )
  }
}

# ------------------------- 19. 画 sample ages -------------------------

for (i in seq_len(numDemes)) {
  if (sampTimes[i] == 0) {
    mtext(
      sampTimes[i],
      side = 1,
      line = -1.2,
      at = i,
      cex = 0.7,
      col = ageCol
    )
  } else {
    text(
      i,
      sampTimes[i],
      labels = sampTimes[i],
      cex = 0.7,
      pos = 1,
      col = ageCol
    )
  }
}

# ------------------------- 20. 画图例 -------------------------

minY.coo <- grconvertY(0, from = "nic", to = "user")

par(xpd = NA)

maxOrder <- ceiling(log10(maxPopSize))
minOrder <- floor(log10(max(minPopSize, 1)))

if (drawLogPopSize) {
  popSizeRadius <- 10 ^ (maxOrder:minOrder)
} else {
  numSteps <- max(length(maxOrder:minOrder) - 1, 1)
  step <- (maxPopSize - minPopSize) / numSteps
  popSizeRadius <- maxPopSize - step * (0:numSteps)
}

winWidth <- numSamples + 2
ypos <- 3 / 4 * minY.coo

text(
  x = -winWidth / 10 * 1.2,
  y = ypos,
  labels = "Pop. \nsizes ",
  cex = 0.8,
  pos = 2
)

for (i in seq_along(popSizeRadius)) {

  curRadius <- interpolRadius(
    popSizeRadius[i],
    minPopSize,
    maxPopSize,
    minRadius,
    maxRadius,
    drawLogPopSize
  )

  xpos <- -winWidth / 12 + (i - 1) * winWidth / 12

  if (curRadius > 0) {
    draw.circle(
      xpos,
      ypos,
      radius = curRadius,
      col = popCol,
      border = popBorderCol
    )
  }

  text(
    xpos,
    ypos - abs(ypos) * 0.1,
    format(popSizeRadius[i], digits = 1, scientific = TRUE),
    cex = 0.7,
    pos = 1,
    col = "black"
  )
}

if (printPDF) {
  dev.off()
  message("PDF 已生成：", pdfFileName)
} else {
  message("绘图完成。若想输出 PDF，请把 printPDF <- TRUE。")
}
