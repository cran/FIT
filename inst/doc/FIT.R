## ---- eval=FALSE---------------------------------------------------------
#  install.packages('FIT')

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("FIT", INSTALL_opts = "--no-multiarch")

## ----chunk-load, eval=TRUE-----------------------------------------------
requireNamespace('FIT')

## ----chunk-train-data, eval=TRUE-----------------------------------------
train.attribute.file <- system.file('extdata', 'train.attribute', package='FIT')
train.weather.file <- system.file('extdata', 'train.weather', package='FIT')
train.expression.file <- system.file('extdata', 'train.expression', package='FIT')

training.attribute  <- FIT::load.attribute(train.attribute.file);
training.weather    <- FIT::load.weather(train.weather.file, 'weather',
                                        c('temperature', 'radiation'))
training.expression <- FIT::load.expression(train.expression.file, 'ex', 
                                           c('Os12g0189300', 'Os02g0724000', 'Os02g0139700', 'Os06g0133200'))

## ----chunk-grid, eval=TRUE-----------------------------------------------
grid.coords <- list(
  env.temperature.threshold = c(10, 15, 20, 25, 30),
  env.temperature.amplitude = c(-100/30, -1/30, 1/30, 100/30),
  env.radiation.threshold = c(1, 10, 20, 30, 40),
  env.radiation.amplitude = c(-100/80, -1/80, 1/80, 100/80),  
  env.temperature.period = c(10, 30, 90, 270, 720, 1440, 1440*3),
  env.radiation.period = c(10, 30, 90, 270, 720, 1440, 1440*3),
  gate.temperature.phase = seq(0, 23*60, 1*60),
  gate.radiation.phase = seq(0, 23*60, 1*60),
  gate.temperature.threshold = cos(pi*seq(4,24,4)/24),
  gate.radiation.threshold = cos(pi*seq(4,24,4)/24),
  gate.temperature.amplitude = c(-5, 5),
  gate.radiation.amplitude = c(-5, 5)
)

## ----chunk-recipe, eval=TRUE---------------------------------------------
recipe <- FIT::make.recipe(c('temperature', 'radiation'), 
                           init = 'gridsearch',
                           optim = c('lm'),
                           fit = 'fit.lasso',
                           init.data = grid.coords,
                           time.step = 10)

## ----chunk-train, eval=TRUE----------------------------------------------
models <- FIT::train(training.expression,
                     training.attribute,
                     training.weather,
                     recipe)

## ---- eval=TRUE----------------------------------------------------------
models <- unlist(models)

## ----chunk-predict, eval=TRUE--------------------------------------------
prediction.attribute.file <- system.file('extdata', 'prediction.attribute', package = 'FIT')
prediction.weather.file <- system.file('extdata', 'prediction.weather', package = 'FIT')

prediction.attribute  <- FIT::load.attribute(prediction.attribute.file);
prediction.weather    <- FIT::load.weather(prediction.weather.file, 'weather',
                                            c('temperature', 'radiation'))
prediction <- FIT::predict(models, prediction.attribute, prediction.weather)

## ----chunk-error, fig.height=2, eval=TRUE--------------------------------
prediction.expression.file <- system.file('extdata', 'prediction.expression', package = 'FIT')

prediction.expression <- FIT::load.expression(prediction.expression.file, 'ex', 
                                            c('Os12g0189300', 'Os02g0724000', 'Os02g0139700', 'Os06g0133200'))
prediction.errors <- FIT::prediction.errors(models,
                                           prediction.expression,
                                           prediction.attribute,
                                           prediction.weather)

## ----chunk-plot, fig.height=4, eval=TRUE---------------------------------
for(i in 1:length(prediction)){
  plot(prediction[[i]], prediction.expression$rawdata[,i], 
      xlab='prediction', ylab='observation')
  title(models[[i]]$gene)
}

## ---- eval=TRUE----------------------------------------------------------
models[[1]]$coefs

## ---- eval=TRUE----------------------------------------------------------
models[[1]]$params

## ----chunk-load2, eval=TRUE----------------------------------------------
genes <- c('Os03g0197000', 'Os01g0892600', 'Os07g0630800', 'Os01g0700100')
training.expression2 <- FIT::load.expression(train.expression.file, 'ex', genes)

## ----chunk-recipe2, eval=TRUE--------------------------------------------
init.params <- rep(list(models[[1]]$params), 4)
names(init.params) <- genes
recipe2 <- FIT::make.recipe(models[[1]]$env, 
                           init = 'manual',
                           optim = c('lm'),
                           fit = 'fit.lasso',
                           init.data = list(
                             params = init.params,
                             response.type = models[[1]]$response.type,
                             input.mean = models[[1]]$input.mean,
                             input.sd = models[[1]]$input.sd
                             ),
                           time.step = 10)

## ----chunk-example2, fig.height=4, eval=TRUE-----------------------------
models2 <- unlist(FIT::train(training.expression2,
                             training.attribute,
                             training.weather,
                             recipe2))
prediction2 <- FIT::predict(models2, prediction.attribute, prediction.weather)

prediction.expression2 <- FIT::load.expression(prediction.expression.file, 'ex', genes)

for(i in 1:length(prediction2)){
  plot(prediction2[[i]], prediction.expression2$rawdata[,i], 
      xlab='prediction', ylab='observation')
  title(models2[[i]]$gene)
}

## ----chunk-weight, eval=TRUE---------------------------------------------
rna.seq.file <- system.file('extdata', 'rna-seq', package='FIT')
weight <- FIT::load.weight(rna.seq.file, 'weights', genes)

## ----chunk-load-rna-seq, eval=TRUE---------------------------------------
training.expression.rnaseq <- FIT::load.expression(rna.seq.file, 'log.cpm', genes)

## ----chunk-recipe-rna-seq, eval=TRUE-------------------------------------
recipe.rnaseq <- FIT::make.recipe(c('temperature', 'radiation'), 
                           init = 'gridsearch',
                           optim = c('lm'),
                           fit = 'fit.lasso',
                           init.data = grid.coords,
                           time.step = 10)

## ----chunk-train-rna-seq, eval=TRUE--------------------------------------
models.rnaseq <- unlist(FIT::train(training.expression.rnaseq,
                             training.attribute,
                             training.weather,
                             recipe.rnaseq,
                             weight))

## ----chunk-predict-rna-seq, eval=TRUE------------------------------------
prediction.rnaseq <- FIT::predict(models.rnaseq, prediction.attribute, prediction.weather)

for(i in 1:length(prediction.rnaseq)){
  plot(prediction.rnaseq[[i]], prediction.expression2$rawdata[,i], 
      xlab='prediction', ylab='observation')
  title(models.rnaseq[[i]]$gene)
}

