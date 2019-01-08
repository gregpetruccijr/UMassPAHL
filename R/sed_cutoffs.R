#' Determines sedentary periods using various methods
#'
#' \code{sed_cutoffs} Will return the 1-second accelerometer data with appended columns corresponding to the various methods indicating sedentary periods
#'
#' @param ag_data_1sec Trimmed hip accelerometer data from trim.ag
#'
#' @param ag_data_1sec_wrist Trimmed wrist accelerometer data from trim.ag
#' @param ag_data_raw Trimmed raw accelerometer data from trim.ag
#' @param subj_label Subject label identifier that will be appended as the first column in the returned data
#' @param cpm100 Apply 100 counts/minute cutpoint to count data (Matthews et al. 2008)
#' @param cpm150 Apply 150 counts/minute cutpoint to count data (Kozey-Keadle et al. 2011)
#' @param cpm200vm Apply 200 counts/minute using vector magnitude cutpoint to count data (Aguilar-Farias et al. 2014)
#' @param soj1x Apply Sojourn 1x algorithm to count data (Lyden et al. 2014)
#' @param soj3x Apply Sojourn 3x algorithm to count data (Lyden et al. 2014)
#' @param jstaud2015 Apply Decision Tree and Random Forest machine learning algorithms to raw wrist data (Staudenmayer et al. 2015)
#' @param sedsphere Apply the Sedentary Sphere to raw wrist data (Rowlands et al. 2016)
#' @param wcpm1853vm Apply 1853 counts/minute using vector magnitude cutpoint to count data (Koster et al. 2016)
#' @param wcpm15s376vm Apply 376 counts/15 seconds using vector magnitude cutpoint to count data (Koster et al. 2016)
#' @param enmo44.8 Apply 44.8 mg/minute using raw wrist data from read.ag or trim.ag (Hildebrand et al. 2017)
#' @param VMcorrG_mod_15s Vector magnitude corrected for gravity cutpoint used by Sedentary Sphere to determine activity periods, default is 489 mg
#'
#' @return Returns the 1-second accelerometer file with columns appended associated with the various methods to determine sedentary periods
#'
#' @examples
#' sed_data = sed_cutoffs(Session_AG_hip_1sec_data, Session_AG_wrist_1sec_data,Session_AG_wrist_raw_data, cpm100 = T)
#'

sed_cutoffs = function(ag_data_1sec, ag_data_1sec_wrist, ag_data_raw, subj_label = participant,
                       cpm100 = F, cpm150 = F, cpm200vm = F, soj1x = F, soj3x = F,
                       sedsphere = F, jstaud2015 = F,
                       wcpm1853vm = F, wcp15s376vm = F,enmo44.8 = F,
                       VMcorrG_mod_15s = 489){

  # Matthews et al. 2008
  if(cpm100 == T){
    ag_data_60sec = ag_epochr(ag_data_1sec, epoch = 60)
    ag_data_60sec$cpm100 = ifelse(ag_data_60sec$Axis1 >=100,0,1)

    expand_ag_data_60sec = ag_data_60sec[rep(seq_len(nrow(ag_data_60sec)), each=60),]

    # handles 1 sec data that doesn't have complete minute
    if(nrow(ag_data_1sec) > nrow(expand_ag_data_60sec)){
      residual = nrow(ag_data_1sec) %% nrow(expand_ag_data_60sec)

      residual_nan = matrix(NaN, nrow = residual, ncol = ncol(expand_ag_data_60sec))
      colnames(residual_nan) = colnames(expand_ag_data_60sec)
      expand_ag_data_60sec = rbind(expand_ag_data_60sec,residual_nan)
    }

    ag_data_1sec = cbind(ag_data_1sec, cpm100 = expand_ag_data_60sec$cpm100)
  }

  # Kozey-Keadle et al. 2011
  if(cpm150 == T){
    ag_data_60sec = ag_epochr(ag_data_1sec, epoch = 60)
    ag_data_60sec$cpm150 = ifelse(ag_data_60sec$Axis1 >=150,0,1)

    expand_ag_data_60sec = ag_data_60sec[rep(seq_len(nrow(ag_data_60sec)), each=60),]

    if(nrow(ag_data_1sec) > nrow(expand_ag_data_60sec)){
      residual = nrow(ag_data_1sec) %% nrow(expand_ag_data_60sec)

      residual_nan = matrix(NaN, nrow = residual, ncol = ncol(expand_ag_data_60sec))
      colnames(residual_nan) = colnames(expand_ag_data_60sec)
      expand_ag_data_60sec = rbind(expand_ag_data_60sec,residual_nan)
    }

    ag_data_1sec = cbind(ag_data_1sec, cpm150 = expand_ag_data_60sec$cpm150)
  }

  # Aguilar-Farias et al. 2014
  if(cpm200vm == T){
    ag_data_60sec = ag_epochr(ag_data_1sec, epoch = 60)
    ag_data_60sec$cpm200vm = ifelse(ag_data_60sec$VM >=200,0,1)

    expand_ag_data_60sec = ag_data_60sec[rep(seq_len(nrow(ag_data_60sec)), each=60),]

    if(nrow(ag_data_1sec) > nrow(expand_ag_data_60sec)){
      residual = nrow(ag_data_1sec) %% nrow(expand_ag_data_60sec)

      residual_nan = matrix(NaN, nrow = residual, ncol = ncol(expand_ag_data_60sec))
      colnames(residual_nan) = colnames(expand_ag_data_60sec)
      expand_ag_data_60sec = rbind(expand_ag_data_60sec,residual_nan)
    }

    ag_data_1sec = cbind(ag_data_1sec, cpm200vm = expand_ag_data_60sec$cpm200vm)
  }

  # Lyden et al. 2014
  if(soj1x == T){
    sojourn.1x.adapted <- function(counts,perc.cut=0.05,perc.cut.2=0.12,perc.cut.3=0.55,too.short=10,sit.cut=90,long.soj=120)
    {
      y <- counts
      # identify sojourns.
      inds <- 1:length(y)

      mmm <- length(y)
      one <- y[-mmm]
      two <- y[-1]

      # transitions from 0 to >0
      trans.up <- (one==0)&(two>0)
      # transitions from >0 to 0
      trans.down <- (one>0)&(two==0)

      trans <- c(0,trans.up+trans.down)
      trans.inds <- (1:mmm)[trans==1]

      # indices where transitions take place
      trans.inds <- c(1,trans.inds,(mmm+1))

      # how long are the sojourns and the zeros
      durations <- trans.inds[-1]-trans.inds[-length(trans.inds)]

      # identify if interval is zeros or >0s (they alternate)
      type <- rep("zeros",length=length(durations))
      if (y[1]==0)
        type <- rep(c("zeros","act"),length=length(durations))
      if (y[1]>0)
        type <- rep(c("act","zeros"),length=length(durations))

      soj.table <- data.frame(type,durations,trans.inds=trans.inds[-length(trans.inds)])

      soj.table$act.type.1 <- "undetermined"
      soj.table$act.type.1[(soj.table$type=="zeros")&(soj.table$durations>sit.cut)] <- "sedentary"
      soj.table$act.type.1[(soj.table$type=="act")&(soj.table$durations>too.short)] <- "activity"



      # combine neighboring undetermineds
      mmm <- dim(soj.table)[1]
      prev.was.undet.inds <-
        (2:mmm)[(soj.table$act.type.1[2:mmm]=="undetermined")&
                  (soj.table$act.type.1[1:(mmm-1)]=="undetermined")]
      if (length(prev.was.undet.inds)>0){
        rev.soj.table <- soj.table[-prev.was.undet.inds,]
      } else {
        rev.soj.table = soj.table
      }

      mmm <- dim(rev.soj.table)[1]

      rev.soj.table$durations <-
        c((rev.soj.table$trans.inds[-1]-
             rev.soj.table$trans.inds[-mmm]),
          rev.soj.table$durations[mmm])

      mmm <- dim(rev.soj.table)[1]

      # find too short undetermineds
      too.short.undet.inds <- (1:mmm)[(rev.soj.table$durations<too.short)&(rev.soj.table$act.type.1=="undetermined")]

      if (length(too.short.undet.inds)>0)
      {
        while (too.short.undet.inds[1]==1)
        {
          too.short.undet.inds <- too.short.undet.inds[-1]
          rev.soj.table <- rev.soj.table[-1,]
          rev.soj.table$trans.inds[1] <- 1
          mmm <- dim(rev.soj.table)[1]
          too.short.undet.inds <- too.short.undet.inds-1
        }

        last <- length(too.short.undet.inds)
        if(last == 1){
          too.short.undet.inds <- too.short.undet.inds[-last]
          junk <- rev.soj.table$durations[(mmm-1)]
          rev.soj.table <- rev.soj.table[-mmm,]
          mmm <- dim(rev.soj.table)[1]
          rev.soj.table$durations[mmm] <- junk+rev.soj.table$durations[mmm]
        } else {
          while (too.short.undet.inds[last]==mmm)
          {
            too.short.undet.inds <- too.short.undet.inds[-last]
            junk <- rev.soj.table$durations[(mmm-1)]
            rev.soj.table <- rev.soj.table[-mmm,]
            mmm <- dim(rev.soj.table)[1]
            rev.soj.table$durations[mmm] <- junk+rev.soj.table$durations[mmm]
            last <- length(too.short.undet.inds)
          }
        }

        # short undetermineds between two acts of same type
        to.delete.inds <-
          (too.short.undet.inds)[rev.soj.table$act.type.1[too.short.undet.inds-1]==rev.soj.table$act.type.1[too.short.undet.inds+1]]
        done.inds <- (1:length(too.short.undet.inds))[rev.soj.table$act.type.1[too.short.undet.inds-1]==rev.soj.table$act.type.1[too.short.undet.inds+1]]
        too.short.undet.inds <- too.short.undet.inds[-done.inds]

        # between two acts of different types
        junk <- rev.soj.table[too.short.undet.inds,]

        if(nrow(junk) > 0){
          junk$act.type.1 <- "sedentary"
          junk$act.type.1[junk$type=="act"] <- "activity"
          rev.soj.table[too.short.undet.inds,] <- junk
        }

        if(length(to.delete.inds) != 0){
          rev.soj.table <- rev.soj.table[-to.delete.inds,]
        }

      }


      mmm <- dim(rev.soj.table)[1]
      junk <- c(rev.soj.table$act.type.1[2:mmm]==rev.soj.table$act.type.1[1:(mmm-1)])
      same.as.prev.inds <- (2:mmm)[junk]
      if (length(same.as.prev.inds)>0)
      {
        rev.soj.table <- rev.soj.table[-same.as.prev.inds,]
        mmm <- dim(rev.soj.table)[1]
        rev.soj.table$durations <-
          c((rev.soj.table$trans.inds[-1]-
               rev.soj.table$trans.inds[-mmm]),
            rev.soj.table$durations[mmm])
        last.obs <- rev.soj.table$durations[mmm]-1+rev.soj.table$trans.inds[mmm]

        if (last.obs != length(y))
          rev.soj.table$durations[mmm] <- length(y)-rev.soj.table$trans.inds[mmm]+1

      }

      trans.inds <- c(rev.soj.table$trans.inds,length(y)+1)
      durations <- trans.inds[-1]-trans.inds[-length(trans.inds)]

      soj.table <- data.frame(durations)

      sojourns <- rep(1:length(soj.table$durations),soj.table$durations)
      perc.gt.0 <- tapply(y>0,sojourns,mean)

      soj.table$perc.gt.0 <- perc.gt.0

      soj.table$revised.type <- "sit.still"
      soj.table$revised.type[soj.table$perc.gt.0>perc.cut.3] <- "activity"
      soj.table$revised.type[(soj.table$perc.gt.0>perc.cut)&(soj.table$perc.gt.0<=perc.cut.2)&(soj.table$durations>sit.cut)] <- "sit.move"
      soj.table$revised.type[(soj.table$perc.gt.0>perc.cut)&(soj.table$perc.gt.0<=perc.cut.2)&(soj.table$durations<=sit.cut)] <- "stand.still"
      soj.table$revised.type[(soj.table$perc.gt.0>perc.cut.2)&(soj.table$perc.gt.0<=perc.cut.3)] <- "stand.small.move"

      durations <- soj.table$durations
      type <- soj.table$revised.type

      sojourns <- rep(1:length(durations),durations)
      type <- rep(type,durations)
      perc.gt.0 <- rep(perc.gt.0,durations)
      durations <- rep(durations,durations)
      nnn <- length(sojourns)

      longer.acts <- unique(sojourns[(durations>(long.soj-1))])

      f <- function(s)
      {
        dur <- 	unique(durations[sojourns==s])
        sub.sojourns <- rep(1:floor(dur/(long.soj/2)),
                            times=c(rep((long.soj/2),floor(dur/(long.soj/2))-1),
                                    dur-(floor(dur/(long.soj/2))-1)*(long.soj/2)))
        sub.sojourns <- s + sub.sojourns/(max(sub.sojourns)+1)
        return(sub.sojourns)
      }
      new.values <- sapply(longer.acts,f)
      starts <- sapply(match(longer.acts,sojourns),paste,":",sep="")
      ends <- length(sojourns) - match(longer.acts,rev(sojourns)) + 1
      indices <- mapply(paste,starts,ends,MoreArgs=list(sep=""),USE.NAMES=F)
      indices <- unlist(lapply(parse(text = indices), eval))
      sojourns[indices] <- unlist(new.values)

      # apply METs to zeros
      METs <- rep(NA,length(type))
      METs[(type=="sit.still")] <- 1
      METs[(type=="sit.move")] <- 1.2
      METs[(type=="stand.still")] <- 1.5
      METs[(type=="stand.small.move")] <- 1.7


      data <- data.frame(counts=y,sojourns=sojourns,durations=durations,type=type,METs=METs,perc.gt.0=perc.gt.0)

      # prepare to apply nnet to the activity sojourns
      nnn <- dim(data)[1]
      act.inds <- (1:nnn)[(data$type=="activity")]
      if(length(act.inds) > 0){
        act.data <- data[act.inds,]
        act.durations <- table(act.data$sojourns)

        quantiles <- tapply(act.data$counts,act.data$sojourns,quantile,p=c(.1,.25,.5,.75,.9))
        nn.data <- as.data.frame(do.call("rbind",quantiles))
        nn.data$acf <- tapply(act.data$counts,act.data$sojourns,acf.lag1)
        nn.data <- nn.data[,c(1:6)]

        names(nn.data) <- c("X10.","X25.","X50.","X75.","X90.","acf")

        nnetinputs <- scale(nn.data,center=cent,scale=scal)

        # apply nnet and put it back into the dataset
        est.mets.1 <- NA #predict(MA.reg.nn,nnetinputs)
        est.mets.2 <- predict(ALL.reg.nn,nnetinputs)

        #act.mets.1 <- rep(est.mets.1,act.durations)
        act.mets.2 <- rep(est.mets.2,act.durations)

        data$METs <- METs
        data$METs.2 <- METs

        data$METs[act.inds] <- act.mets.2
        data$METs.2[act.inds] <- act.mets.2
      }
      data$level <- "sed"
      data$level[data$METs>=1.5] <- "light"
      data$level[data$METs>=3] <- "mod"
      data$level[data$METs>=6] <- "vig"
      data$level <- factor(data$level,levels=c("sed","light","mod","vig"))

      data$level.2 <- "sed"
      data$level.2[data$METs.2>=1.5] <- "light"
      data$level.2[data$METs.2>=3] <- "mod"
      data$level.2[data$METs.2>=6] <- "vig"
      data$level.2 <- factor(data$level.2,levels=c("sed","light","mod","vig"))
      n <- dim(data)[1]
      inds <- (1:n)[data$METs<1]
      data$METs[inds] <- 1

      data <- data[,c(1,2,3,4,5,6,8)]
      return(data)
    }

    soj1x_data = sojourn.1x.adapted(ag_data_1sec$Axis1)

    soj1x_data$soj1x = ifelse(soj1x_data$METs < 1.5, 1,0)

    #soj1x_sed_indices = which(str_detect(soj1x_data$type, 'sit'))

    ## Below is an sed/not designation based on MET calculations. Check w/ JStaud if MET designation or thresholds would be best
    #soj1x_data$soj1x = ifelse(soj1x_data[,7] == 'sed', 'sedentary','non-sedentary')

    ag_data_1sec = cbind(ag_data_1sec, soj1x = soj1x_data$soj1x)

  }

  if(soj3x == T){
    sojourn.3x.adapted <- function(counts,counts.2,counts.3,vect.mag,short=30)
    {
      y <- counts
      counts.2 <- counts.2
      counts.3 <- counts.3

      inds <- 1:length(y)

      mmm <- length(y)
      one <- y[-mmm]
      two <- y[-1]

      # find transitions

      ###### Break due to no transitions
      trans <- ((one-two)>15)&(two<=10) 	# this is how i find initial transitions

      if(sum(trans) > 0){
        trans <- c(0,trans)

        trans.inds <- (1:mmm)[trans==1]

        # how long between transistions

        durations <- trans.inds[-1]-trans.inds[-length(trans.inds)]

        #	put first duration in and make last trans go till end of file

        dd <- length(durations)
        tt <- length(trans.inds)
        durations[dd+1] <- mmm-trans.inds[tt]
        dd <- length(durations)
        durations.junk <- trans.inds[1]
        durations <- c(durations.junk,durations)
        dd <- length(durations)

        durations.compare <- durations
        length(durations.compare)

        # get number of sojourns

        sojourns <- rep(1:length(durations),1)
        sojourns.long <- rep(sojourns,durations)
        mean.cnts.soj <- as.vector(tapply(y,sojourns.long,mean))

        # combine too short sojourns.

        #	combine too short sojourns with neighboring sojourn.
        # 	this loop repeats until there are no more too short sojourns

        counter <- 1

        repeat	# loop 1

        {
          too.short <- (1:dd)[durations<short]
          ts <- length(too.short)

          if(length(too.short)==0)
            break

          if(length(too.short)>0)
          {


            # this loop deals with instances where the first too.short sojourn is first sojourn of file ie. it only has a second neighbor to combine it with

            counter.1 <- 1

            repeat	 # loop 2

            {

              if (too.short[counter.1]==counter.1)
              {
                sojourns[1:counter.1] <- sojourns[counter.1+1]

                counter.1 <- counter.1+1
              }

              if(is.na(too.short[counter.1])){
                break
              } else {
                if (too.short[counter.1]!=counter.1)

                  break
              }
            }	# end loop 2

            s <- length(sojourns)

            # this loop deals with if last too short sojourn is last sojourn of file ie. it only has a first neighbor to combine it with

            counter.2 <- s
            counter.ts <- ts

            repeat{

              if (too.short[counter.ts]==counter.2)
              {
                sojourns[counter.2:s] <- sojourns[counter.2-1]

                counter.2 <- counter.2-1
                counter.ts <- counter.ts-1
              }

              if (too.short[counter.ts]!=counter.2)

                break

            }	#end loop 3

            s <- length(sojourns)

            # now deal with all other too short sojourns

            junk.too.short <- too.short

            if(counter.ts<ts-1)
            {
              junk.too.short <- too.short[-(counter.ts+1:ts)]
            }
            if (counter.1>1)
            {
              junk.too.short <- junk.too.short[-(1:counter.1-1)]
            }

            j.t.s <- length(junk.too.short)

            first.neighbors <- junk.too.short-1
            second.neighbors <- junk.too.short+1

            #	right now i combine too short sojourns with its neighbor that was shorter in duration (e.g. first neighbor = 60 seconds long and second neighbor = 300 seconds long, it gets combined with first neighbor)

            revised.sojourns <- sojourns

            durations[junk.too.short]

            durations.first.neighbors <- durations[first.neighbors]
            durations.second.neighbors <- durations[second.neighbors]

            #	put in dummy duration for too.short sojourns at beginning and end of file
            durations.first.neighbors[is.na(durations.first.neighbors)] <- 100000
            durations.second.neighbors[is.na(durations.second.neighbors)] <- 100000

            n.neighbors <- length(durations.first.neighbors)
            n.neighbors.2 <- length(durations.second.neighbors)

            inds.first <- (1:n.neighbors)[durations.first.neighbors<=durations.second.neighbors]
            inds.second <- (1:n.neighbors)[durations.first.neighbors>durations.second.neighbors]

            too.short.inds.first <- junk.too.short[inds.first]
            too.short.inds.second <- junk.too.short[inds.second]

            revised.sojourns[too.short.inds.first] <- first.neighbors[inds.first]
            revised.sojourns[too.short.inds.second] <- second.neighbors[inds.second]

            # deal with instances where need to combine more than 2 sojourns - i.e. short sojourn became first neighbor, and then sojourn before first neighbor also becomes that sojourn via second neighbor grouping - want all 3 of these sojourns to be combined.

            rs <- length(revised.sojourns)

            one.order <- revised.sojourns[-rs]
            two.order <- revised.sojourns[-1]

            o <- length(one.order)

            inds.order <- (1:o)[one.order>two.order]
            if (length(inds.order>0))
              revised.sojourns[inds.order+1] <- revised.sojourns[inds.order]

            # get new durations now that sojourns are combined

            rs <- length(revised.sojourns)
            revised.durations <- as.vector(tapply(durations,revised.sojourns,sum))

            rd <- length(revised.durations)

            # get new sojourns now that durations are combined

            revised.sojourns <- rep(1:length(revised.durations),1)
            rs <- length(revised.sojourns)

            durations <- revised.durations
            dd <- length(durations)
            sojourns <- revised.sojourns
            s <- length(sojourns)

          }

          #	print(counter)
          counter <- counter+1

        }	# end loop 1

      } else {
        durations = length(y)
        sojourns = 1
      }
      #	 make table of durations and sojourns etc

      trans.table <- data.frame(counts=y,counts.2=counts.2,counts.3=counts.3,vect.mag=vect.mag,sojourns=0,durations=0,perc.soj=NA,soj.type.all=NA,soj.mets.all=NA)

      tt <- dim(trans.table)[1]
      durations.1 <- rep(durations,durations)
      sojourns.1 <- rep(sojourns,durations)

      trans.table$durations <- durations.1
      trans.table$sojourns <- sojourns.1

      #	get percent non zero in table

      perc.soj <- tapply(y>0,sojourns.1,mean)

      perc.soj <- rep(perc.soj,durations)

      trans.table$perc.soj <- perc.soj


      ### get inds.inactivities so can test nnet only to distinguish between lifestyle and sedentary

      #	now get inactivity indices

      inds.inacts <- (1:tt)[trans.table$perc.soj<0.7]
      inactivities <- trans.table[inds.inacts,]
      i.a <- dim(inactivities)[1]

      inact.trans.inds <- c(1,(1+(1:i.a)[inactivities$sojourns[-1]!=inactivities$sojourns[-i.a]]))

      inact.durations <- inactivities$durations[inact.trans.inds]

      #	get nnetinputs for vertical axis

      nnetinputs <-
        as.vector(unlist(tapply(inactivities$counts,inactivities$sojourns,quantile,probs=c(.1,.25,.5,.75,.9))))
      length(nnetinputs)
      nnetinputs <- matrix(nnetinputs,length(nnetinputs)/5,5,byrow=T)
      nnetinputs <- as.data.frame(nnetinputs)
      names(nnetinputs) <- c("X10.","X25.","X50.","X75.","X90.")
      nnetinputs$acf <- 0

      g <- 1
      for (soj in unique(inactivities$sojourns))
      {
        counts <- inactivities$counts[inactivities$sojourns==soj]


        if (sum(counts)>0)
        {
          temp <- acf(counts,lag.max=1,plot=F)
          nnetinputs$acf[g] <- as.numeric(unlist(temp[1,1])[1])

        }
        g <- g+1
        #	print(g)
      }

      nnetinputs$acf[is.na(nnetinputs$acf)] <-
        mean(nnetinputs$acf,na.rm=T)

      ####	get nnetinputs.2 - second axis

      nnetinputs.2 <-
        as.vector(unlist(tapply(inactivities$counts.2,inactivities$sojourns,quantile,probs=c(.1,.25,.5,.75,.9))))
      length(nnetinputs.2)
      nnetinputs.2 <- matrix(nnetinputs.2,length(nnetinputs.2)/5,5,byrow=T)
      nnetinputs.2 <- as.data.frame(nnetinputs.2)
      names(nnetinputs.2) <- c("X10.2","X25.2","X50.2","X75.2","X90.2")
      nnetinputs.2$acf.2 <- 0

      g <- 1
      for (soj in unique(inactivities$sojourns))
      {
        counts <- inactivities$counts.2[inactivities$sojourns==soj]


        if (sum(counts)>0)
        {
          temp <- acf(counts,lag.max=1,plot=F)
          nnetinputs.2$acf.2[g] <- as.numeric(unlist(temp[1,1])[1])

        }
        g <- g+1
        #	print(g)
      }

      nnetinputs.2$acf.2[is.na(nnetinputs.2$acf.2)] <-
        mean(nnetinputs.2$acf.2,na.rm=T)


      ####get nnetinputs.3 - third axis

      nnetinputs.3 <-
        as.vector(unlist(tapply(inactivities$counts.3,inactivities$sojourns,quantile,probs=c(.1,.25,.5,.75,.9))))
      length(nnetinputs.3)
      nnetinputs.3 <- matrix(nnetinputs.3,length(nnetinputs.3)/5,5,byrow=T)
      nnetinputs.3 <- as.data.frame(nnetinputs.3)
      names(nnetinputs.3) <- c("X10.3","X25.3","X50.3","X75.3","X90.3")
      nnetinputs.3$acf.3 <- 0

      g <- 1
      for (soj in unique(inactivities$sojourns))
      {
        counts <- inactivities$counts.3[inactivities$sojourns==soj]


        if (sum(counts)>0)
        {
          temp <- acf(counts,lag.max=1,plot=F)
          nnetinputs.3$acf.3[g] <- as.numeric(unlist(temp[1,1])[1])

        }
        g <- g+1
        #print(g)
      }

      nnetinputs.3$acf.3[is.na(nnetinputs.3$acf.3)] <-
        mean(nnetinputs.3$acf.3,na.rm=T)

      ####get nnetinputs.vm - vector magnitude

      nnetinputs.vm <-
        as.vector(unlist(tapply(inactivities$vect.mag,inactivities$sojourns,quantile,probs=c(.1,.25,.5,.75,.9))))
      length(nnetinputs.vm)
      nnetinputs.vm <- matrix(nnetinputs.vm,length(nnetinputs.vm)/5,5,byrow=T)
      nnetinputs.vm <- as.data.frame(nnetinputs.vm)
      names(nnetinputs.vm) <- c("X10.vm","X25.vm","X50.vm","X75.vm","X90.vm")
      nnetinputs.vm$acf.vm <- 0

      g <- 1
      for (soj in unique(inactivities$sojourns))
      {
        counts <- inactivities$vect.mag[inactivities$sojourns==soj]


        if (sum(counts)>0)
        {
          temp <- acf(counts,lag.max=1,plot=F)
          nnetinputs.vm$acf.vm[g] <- as.numeric(unlist(temp[1,1])[1])

        }
        g <- g+1
        #print(g)
      }

      nnetinputs.vm$acf.vm[is.na(nnetinputs.vm$acf.vm)] <-
        mean(nnetinputs.vm$acf.vm,na.rm=T)

      #	combine inputs so can center and scale

      inputs <- cbind(nnetinputs,nnetinputs.2)
      inputs <- cbind(inputs,nnetinputs.3)
      inputs <- cbind(inputs,nnetinputs.vm)
      inputs <- cbind(inputs,inact.durations)

      inputs <- scale(inputs,center=cent.1,scale=scal.1)
      inputs <- as.data.frame(inputs)

      #	predict type using all axes + vm.  i intially had a lot of prediction nnets here (ie different axis) but have removed them and only include the one that looks "the best".  there are definitely others we can use/try

      #	remove NA's

      inputs.1 <- inputs[,-(13)]
      inputs.1 <- inputs.1[,-(1:2)]

      cool.all <- predict(class.nnn.6,inputs.1)

      #	add soj.type to trans table

      junk.cool.all <- as.vector(apply(cool.all,1,which.max))

      cool.all <- rep(junk.cool.all,inact.durations)

      trans.table$soj.type.all[inds.inacts] <- cool.all
      #	assign mets to types.

      trans.table$soj.mets.all[(trans.table$soj.type.all==1)&(trans.table$perc.soj<=0.12)] <- 1.5
      trans.table$soj.mets.all[(trans.table$soj.type.all==1)&(trans.table$perc.soj>0.12)] <- 1.7
      trans.table$soj.mets.all[(trans.table$soj.type.all==3)&(trans.table$perc.soj<=0.05)] <- 1
      trans.table$soj.mets.all[(trans.table$soj.type.all==3)&(trans.table$perc.soj>0.05)] <- 1.2

      #	this identifies activities for nnet all - 6 means activity

      trans.table$soj.type.all[trans.table$perc.soj>=0.7] <- 6

      inds.activity.all <- (1:tt)[(trans.table$perc.soj>=0.7)|(trans.table$soj.type.all==2)|(trans.table$soj.type.all==4)]

      act.trans.table.all <- trans.table[inds.activity.all,]
      dim(act.trans.table.all)
      activity.durations.all <- table(act.trans.table.all$sojourns)

      quantiles.all <- tapply(act.trans.table.all$counts,act.trans.table.all$sojourns,quantile,p=c(.1,.25,.5,.75,.9))
      if(length(quantiles.all)>0){
        nn.trans.table.all <- as.data.frame(do.call("rbind",quantiles.all))

        #	i realize i am getting lag1 differently than i do for inactivities...i should change to use function throughout.
        nn.trans.table.all$acf <- tapply(act.trans.table.all$counts,act.trans.table.all$sojourns,acf.lag1)
        nn.trans.table.all <- nn.trans.table.all[,c(1:6)]

        names(nn.trans.table.all) <- c("X10.","X25.","X50.","X75.","X90.","acf")

        nnetinputs.acts.all <- scale(nn.trans.table.all,center=cent,scale=scal)

        #	predict METs

        act.mets.all <- predict(reg.nn,nnetinputs.acts.all)
        act.mets.all <- rep(act.mets.all,activity.durations.all)

        #	put back in table

        trans.table$soj.mets.all[inds.activity.all] <- act.mets.all
      }
      #	get breaks from sitting

      #	trans.table$do.breaks <- 0
      trans.table$soj.breaks.all <- 0


      soj.posture <- as.vector(trans.table$soj.mets.all)
      s.p <- length(soj.posture)

      soj.one.posture <- soj.posture[-s.p]
      soj.two.posture <- soj.posture[-1]

      soj.trans <- (soj.one.posture<1.5)&(soj.two.posture>=1.5)
      soj.trans <- c(0,soj.trans)
      soj.trans.inds <- (1:s.p)[soj.trans==1]

      trans.table$soj.breaks.all <- soj.trans
      #	sum(trans.table$soj.breaks.all)


      names(trans.table)[8:10] <- c("type","METs","break")

      trans.table <- trans.table[,-c(8,10)]

    }	#	end sojourn

    soj3x_data = sojourn.3x.adapted(ag_data_1sec$Axis1, ag_data_1sec$Axis2, ag_data_1sec$Axis3, ag_data_1sec$VM)

    soj3x_data$soj3x = ifelse(soj3x_data$METs < 1.5, 1,0)

    ag_data_1sec = cbind(ag_data_1sec, soj3x = soj3x_data$soj3x)

  }

  # Staudenmayer et al. 2015
  if(jstaud2015 == T){
    win.width <- 15

    n <- dim(ag_data_raw)[1]

    mins <- ceiling(n/(80*win.width))

    ag_data_raw$min <- rep(1:mins,each=win.width*80)[1:n]
    ag_data_raw$v.ang <- 90*(asin(ag_data_raw$AxisX/ag_data_raw$VM)/(pi/2))
    ag_data_raw.sum <- data.frame(mean.vm=tapply(ag_data_raw$VM,ag_data_raw$min,mean,na.rm=T),
                                  sd.vm=tapply(ag_data_raw$VM,ag_data_raw$min,sd,na.rm=T),
                                  mean.ang=tapply(ag_data_raw$v.ang,ag_data_raw$min,mean,na.rm=T),
                                  sd.ang=tapply(ag_data_raw$v.ang,ag_data_raw$min,sd,na.rm=T),
                                  p625=tapply(ag_data_raw$VM,ag_data_raw$min,pow.625),
                                  dfreq=tapply(ag_data_raw$VM,ag_data_raw$min,dom.freq),
                                  ratio.df=tapply(ag_data_raw$VM,ag_data_raw$min,frac.pow.dom.freq))

    # sedentary or not estimates (rf and tree)
    ag_data_raw.sum$sed.rf <- predict(rf.sed.model,newdata=ag_data_raw.sum)
    ag_data_raw.sum$sed.tr <- predict(tr.sed.model,newdata=ag_data_raw.sum,type="class")

    sed.rf = rep(ag_data_raw.sum$sed.rf, each = win.width)
    sed.tr = rep(ag_data_raw.sum$sed.tr, each = win.width)

    ag_data_1sec = cbind(ag_data_1sec, jstaud2015_sed.tr = sed.tr, jstaud2015_sed.rf = sed.rf)

  }

  # Rowlands et al. 2016
  if(sedsphere == T){
    win.width <- 15

    n <- dim(ag_data_raw)[1]

    mins <- ceiling(n/(80*win.width))

    ag_data_raw$min <- rep(1:mins,each=win.width*80)[1:n]

    ag_data_raw.sum <- data.frame(mean.x=tapply(ag_data_raw$AxisX,ag_data_raw$min,mean,na.rm=T),
                                  mean.y=tapply(ag_data_raw$AxisY,ag_data_raw$min,mean,na.rm=T),
                                  mean.z=tapply(ag_data_raw$AxisZ,ag_data_raw$min,mean,na.rm=T),
                                  sum.VMcorrG = tapply(ag_data_raw$VMcorrG,ag_data_raw$min,sum,na.rm=T),
                                  sum.vm = tapply(ag_data_raw$VM, ag_data_raw$min, sum, na.rm = T))

    ag_data_raw.sum$v.ang <- ifelse(ag_data_raw.sum$mean.y > 1, asin(1)*180/pi,
                                    ifelse(ag_data_raw.sum$mean.y < -1, asin(-1)*180/pi,
                                           asin(ag_data_raw.sum$mean.y)*180/pi))

    # 0 = Sedentary, 1 = Standing, 2 = Activity
    ag_data_raw.sum$sedsphere = ifelse(ag_data_raw.sum$sum.VMcorrG > VMcorrG_mod_15s,2,
                                       ifelse(ag_data_raw.sum$v.ang < -15,1,0))

    sedsphere = rep(ag_data_raw.sum$sedsphere, each = win.width)

    ag_data_1sec = cbind(ag_data_1sec, sedsphere = sedsphere)

  }

  # Koster et al. 2016
  if(wcpm1853vm == T){
    ag_data_60sec = ag_epochr(ag_data_1sec_wrist, epoch = 60)
    ag_data_60sec$wcpm1853vm = ifelse(ag_data_60sec$VM >=1853,0,1)

    expand_ag_data_60sec = ag_data_60sec[rep(seq_len(nrow(ag_data_60sec)), each=60),]

    if(nrow(ag_data_1sec_wrist) > nrow(expand_ag_data_60sec)){
      residual = nrow(ag_data_1sec_wrist) %% nrow(expand_ag_data_60sec)

      residual_nan = matrix(NaN, nrow = residual, ncol = ncol(expand_ag_data_60sec))
      colnames(residual_nan) = colnames(expand_ag_data_60sec)
      expand_ag_data_60sec = rbind(expand_ag_data_60sec,residual_nan)
    }

    ag_data_1sec = cbind(ag_data_1sec, wcpm1853vm = expand_ag_data_60sec$wcpm1853vm)
  }

  if(wcp15s376vm == T){
    ag_data_15sec = ag_epochr(ag_data_1sec_wrist, epoch = 15)
    ag_data_15sec$wcp15s376vm = ifelse(ag_data_15sec$VM >=376,0,1)

    expand_ag_data_15sec = ag_data_15sec[rep(seq_len(nrow(ag_data_15sec)), each=15),]

    if(nrow(ag_data_1sec_wrist) > nrow(expand_ag_data_15sec)){
      residual = nrow(ag_data_1sec_wrist) %% nrow(expand_ag_data_15sec)

      residual_nan = matrix(NaN, nrow = residual, ncol = ncol(expand_ag_data_15sec))
      colnames(residual_nan) = colnames(expand_ag_data_15sec)
      expand_ag_data_15sec = rbind(expand_ag_data_15sec,residual_nan)
    }

    ag_data_1sec = cbind(ag_data_1sec, wcp15s376vm = expand_ag_data_15sec$wcp15s376vm)
  }

  # Hildebrand et al. 2017
  if(enmo44.8 == T){
    win.short = 1
    win.width = 60

    n <- dim(ag_data_raw)[1]

    mins.short = ceiling(n/(80*win.short))
    mins <- ceiling(n/(80*win.width))

    ag_data_raw$min.short <- rep(1:mins.short,each=win.short*80)[1:n]
    ag_data_raw$min <- rep(1:mins,each=win.width*80)[1:n]

    ag_data_raw.enmo_1sec = data.frame(mean.enmo=tapply(ag_data_raw$ENMO,ag_data_raw$min.short,mean,na.rm=T),
                                       min = tapply(ag_data_raw$min,ag_data_raw$min.short,mean,na.rm=T))


    ag_data_raw.sum <- data.frame(mean.enmo = tapply(ag_data_raw.enmo_1sec$mean.enmo,ag_data_raw.enmo_1sec$min,mean,na.rm=T))

    ag_data_raw.sum = ag_data_raw.sum*10^3

    ag_data_raw.sum$enmo44.8 = ifelse(ag_data_raw.sum$mean.enmo >=44.8,0,1)

    enmo44.8 = rep(ag_data_raw.sum$enmo44.8, each = win.width)

    if(nrow(ag_data_1sec) > length(enmo44.8)){
      residual = nrow(ag_data_raw) %% nrow(enmo44.8)

      residual_nan = matrix(NaN, nrow = residual, ncol = ncol(enmo44.8))
      colnames(residual_nan) = colnames(enmo44.8)
      enmo44.8 = rbind(enmo44.8,residual_nan)
    } else {
      if(nrow(ag_data_1sec) < length(enmo44.8)){
        enmo44.8 = enmo44.8[1:nrow(ag_data_1sec)]
      }
    }

    ag_data_1sec = cbind(ag_data_1sec, enmo44.8 = enmo44.8)

  }

  ag_data_1sec = cbind(Participant = rep(as.character(subj_label), nrow(ag_data_1sec)), ag_data_1sec, stringsAsFactors= F)

  return(ag_data_1sec)
}

