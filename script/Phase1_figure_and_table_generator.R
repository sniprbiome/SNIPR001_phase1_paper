if(!c("ggplot2") %in% installed.packages()) {install.packages("ggplot2")}
if(!c("cowplot") %in% installed.packages()) {install.packages("cowplot")}
if(!c("reshape2") %in% installed.packages()) {install.packages("reshape2")}
if(!c("RColorBrewer") %in% installed.packages()) {install.packages("RColorBrewer")}
if(!c("dplyr") %in% installed.packages()) {install.packages("dplyr")}
if(!c("tidyr") %in% installed.packages()) {install.packages("tidyr")}
if(!c("openxlsx") %in% installed.packages()) {install.packages("openxlsx")}
if(!c("ggnewscale") %in% installed.packages()) {install.packages("ggnewscale")}


library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggnewscale)



#a function to reverse the output of table to the input for table
untable <- function(vector){
  listed <- lapply(X = seq(1,length(vector)),FUN = function(x,vector){
    return(rep(names(vector)[x],vector[x]))
  },"vector"=vector)
  return(unlist(listed))
}

#calculates the geometric mean
geomean <- function(x){
  exp(mean(log(x)))
}

#A function to easily set axis text to scientific notation
to_tenth_power_labels <- function(x_labels){
  return(parse(text=paste0("10^",x_labels)))
}

#let's try to figure out a way to more intelligently dodge pointranges that are sequential
intelligently_dodge_one_timepoint <- function(point_data,close_dodge,whisker_point_dodge,whisker_dodge){
  
  colnames(point_data) <- c('group','mean','upper','lower')
  
  #we generate every combination worth testing
  combs <- t(combn(levels(point_data[["group"]]),2))
  colnames(combs) <- c('point1','point2')
  
  #we then figure out how close the means are to each other. If the points are
  #very close to each other then there is a point overlap and we will need to offset
  #things more.
  dists_matrix <- as.matrix(dist(point_data[["mean"]],diag = T,upper = T))
  colnames(dists_matrix) <- levels(point_data[["group"]])
  rownames(dists_matrix) <- levels(point_data[["group"]])
  
  #we get if there are overlaps in the whiskers
  whisker_overlaps <- apply(combs,1,FUN = function(x,point_data){
    #we get the relevant data
    pointA_data <- point_data[point_data[["group"]]== x[[1]],c("upper","lower")]
    pointB_data <- point_data[point_data[["group"]]==x[[2]],c("upper","lower")]
    
    #and we check if the whiskers overlap
    pointA_data[["lower"]]<=pointB_data["upper"] & pointA_data[["upper"]]>=pointB_data["lower"]
    
    
  },'point_data'=point_data)
  
  #there is actually a third way of overlap - when the whisker intersects with the point and just the point
  #so we'll check for that
  whisker_point_overlap <- apply(combs,1,FUN = function(x,point_data){
    #we get the relevant data
    pointA_data <- point_data[point_data[["group"]]== x[[1]],c("upper","lower","mean")]
    pointB_data <- point_data[point_data[["group"]]==x[[2]],c("upper","lower","mean")]
    
    #and we check if the whiskers overlap
    pointA_data[["mean"]]>=pointB_data[["lower"]]&pointA_data[["mean"]]<=pointB_data[["upper"]]|
      pointB_data[["mean"]]>=pointA_data[["lower"]]&pointB_data[["mean"]]<=pointA_data[["upper"]]
    
    
  },'point_data'=point_data)
  
  
  #we add together the data to see if there is overlap in whiskers and/or points
  distances <- cbind(as.data.frame(combs),
                     'whisker_overlap'=whisker_overlaps,
                     'whisker_point_overlap'=whisker_point_overlap,
                     'point_overlap'=as.numeric(dists_matrix[combs])<=close_dodge)
  
  #since the way combn (thankfully) does things 1->2, 1->3, 1->4, 2->3, etc we can 
  #set 1 to be "stationary" and simply move points relative to points 1
  #so the way we will solve this is that we will move y and subsequent points
  #if x and y has a conflict and y is the latter of the two in the sequence we use
  #so we will simply compare things one to one until the distance required between two points
  #has been achieved
  offsets <- rep(0,length(levels(point_data[["group"]])))
  names(offsets) <- levels(point_data[["group"]])
  
  for(row in rownames(distances)){
    #we get the overlap type. If the value is 
    #1 we know we have a whisker overlap if it's 
    #2 then we know we have a point-whisker overlap if it's
    #3 then we have a point-point overlap 
    overlap_type <- rowSums(distances[row,c("whisker_overlap","point_overlap","whisker_point_overlap")])
    
    #we get the points as handy dandy variables
    point1 <- distances[row,"point1"]
    point2 <- distances[row,"point2"]
    
    #we get the current distnance between the two points
    current_overlap <- offsets[point2]-offsets[point1]
    
    #then we calculate the overlap we need
    if(overlap_type==1){
      #we add the required overlap between the intended and what we have
      required_overlap=whisker_dodge-current_overlap
    }else if(overlap_type==2){
      required_overlap=whisker_point_dodge-current_overlap
      
    }else if(overlap_type==3){
      #we add the required overlap between the intended and what we have
      required_overlap=close_dodge-current_overlap
      
    }else{
      #if there is no overlap then why even dodge we just stack!
      required_overlap=0
    }
    
    #we don't want something that has been previously dodged to be un-dodged because there are two whiskers next to each other
    if(required_overlap<0){
      required_overlap=0
    }
    
    #then we add that offset to point2 and all the subsequent points
    offsets[seq(which(point2==names(offsets)),length(offsets))]=offsets[seq(which(point2==names(offsets)),length(offsets))]+required_overlap
  }
  
  #we'll center on the range of points
  offsets <- offsets-mean(range(offsets))
  
  return(offsets)
}

#since we'll have to make this type of figure multiple times we'll make a
#function for making the typical boxplot with colors outlier plots and
#by-group FDR-corrected P-values displayed for each day
#col names of the diversity plotdata should be Day, Group, and Diversity
diversity_figure_generator <- function(diversity_plotdata,y_axis_text,treatment_days,cohort_pointshapes,cohort_colors,cohort_scale_labels,custom_theme,treatment_window_highlight_color,highlight_alpha,y_breaks_scale,comparison_type){
  #By default when making boxplots ggplot colors the outlier plots as black and (as far as I know) there is no way to independenty let these
  #outliers vary with a factor without also coloring the border of the boxplot
  #so we'll create data that only includes data outside IQR as defined by Tukey
  dotdat <- diversity_plotdata %>% 
    group_by(Group,Day) %>% 
    mutate(Q1=quantile(Diversity,0.25,na.rm = T),
           Q3=quantile(Diversity,0.75,na.rm = T),
           IQR=Q3-Q1,
           lower.bound=Q1-IQR*1.5,
           upper.bound=Q3+IQR*1.5,
           'outlier'=Diversity>upper.bound|Diversity<lower.bound) %>% 
    filter(outlier)
  
  #and factor the coloring variable
  dotdat[["Group"]] <- factor(dotdat[["Group"]])
  
  #we set the width of the boxplots so we can offset the points appropriately
  boxplot_width <- 0.5
  
  #we make the offset ranges to reflect the mid-point of the four cohorts
  offset_range <- boxplot_width/length(levels(dotdat[["Group"]]))
  group_offset <- seq(offset_range*-(length(unique(dotdat[["Group"]]))/2-0.5),offset_range*(length(unique(dotdat[["Group"]]))/2-0.5),offset_range)
  names(group_offset) <- levels(dotdat[["Group"]])
  
  #then we'll add the x-coordinate to actually the value at as day_as_number
  dotdat <- cbind(dotdat,'day_as_number'=as.numeric(dotdat[["Day"]])+group_offset[dotdat[["Group"]]])
  
  #we'll carry out a bunch of mann-whitney U tests
  diversity_comparisons <- data.frame()
  temp_data_add <- data.frame()
  
  #we define groups to compare the control group to
  noncontrol_groups <- unique(diversity_plotdata[["Group"]])
  noncontrol_groups <- noncontrol_groups[noncontrol_groups!='Placebo']
  
  #We go through each day of comparison
  for(day in levels(diversity_plotdata[["Day"]])){
    
    #and we get the relevant data
    relevant_day_data <- diversity_plotdata[diversity_plotdata[["Day"]]==day,]
    
    #we get the beta diversities for the particular day to compare the other groups with
    control_beta_div <- relevant_day_data[["Diversity"]][relevant_day_data[["Group"]]=='Placebo']
    
    #we go through all the non-control groups
    for(group in noncontrol_groups){
      #and we then get the beta diversity for the non-control group
      group_beta_div <- relevant_day_data[["Diversity"]][relevant_day_data[["Group"]]==group]
      
      #we perform the wilcox test
      test_results <- wilcox.test(control_beta_div,group_beta_div,alternative = comparison_type)
      
      #and save the relevant data
      temp_data_add <- rbind(temp_data_add,data.frame('pval'=test_results[["p.value"]],'group1'='Placebo','group2'=group,'day'=day))
      
    }
    
    #we merge the datasets
    diversity_comparisons <- rbind(diversity_comparisons,temp_data_add)
    temp_data_add <- data.frame()
  }
  
  #we FDR correct on each group since we're testing the same group effect multiple times (day)
  diversity_comparisons <- diversity_comparisons %>% group_by(group2) %>% mutate('adjusted_p'=p.adjust(pval,method = 'fdr'))
  
  #we'll also convert the groups to numerics
  #we sort to make the numerics stuff easier
  diversity_comparisons[["day"]] <- factor(diversity_comparisons[["day"]],levels = as.character(sort(as.numeric(unique(diversity_comparisons[["day"]])))))
  diversity_comparisons <- cbind(diversity_comparisons,
                                 'group1_as_number'=as.numeric(diversity_comparisons[["day"]])+group_offset[diversity_comparisons[["group1"]]],
                                 'group2_as_number'=as.numeric(diversity_comparisons[["day"]])+group_offset[diversity_comparisons[["group2"]]]
  )
  
  
  #we expand the range a bit to ensure we have time for all the P-values and
  #text on top of the actual data
  plot_top <- max(diversity_plotdata[["Diversity"]])+diff(range(diversity_plotdata[["Diversity"]]))*0.4
  plot_bottom <- min(diversity_plotdata[["Diversity"]])-diff(range(diversity_plotdata[["Diversity"]]))*0.1
  
  
  y_axis_range <- plot_top-plot_bottom
  #we'll set the significance data to be on top of each other by a certain increment
  offset_increment <- y_axis_range/20
  
  
  #we'll also set y-values for the comparisons in a stacked manner
  y_heights_offset_vect=c('Cohort 1'=offset_increment,'Cohort 2'=offset_increment*2,'Cohort 3'=offset_increment*3)
  max_vals <- diversity_plotdata %>% group_by(Day) %>% summarise('max_height'=max(Diversity)) %>% as.data.frame
  rownames(max_vals) <- max_vals[["Day"]]
  
  #we add the y axis coordinate
  diversity_comparisons <- cbind(diversity_comparisons,'y'=max_vals[as.character(diversity_comparisons[["day"]]),"max_height"]+y_heights_offset_vect[diversity_comparisons[["group2"]]])
  
  #we'll convert the boxplot data to numeric
  x_breaks <- seq(1,length(levels(diversity_plotdata[["Day"]])))
  x_labels <- levels(diversity_plotdata[["Day"]])
  
  #we get the numbers that are used in treatment days
  treatment_days_numeric <- x_breaks[x_labels%in%treatment_days]
  #and we expand the x-axis labels a bit to each side to include the boxes
  treatment_days_numeric[1] <- treatment_days_numeric[1]-(boxplot_width/2+0.05)
  treatment_days_numeric[length(treatment_days_numeric)] <- treatment_days_numeric[length(treatment_days_numeric)]+(boxplot_width/2+0.05)
  
  #we generate the data to get the 
  ribbon_dat <- data.frame('treatment_days'=treatment_days_numeric,
                           'ymax'=rep(plot_top,length(treatment_days_numeric)),
                           'ymin'=rep(plot_bottom,length(treatment_days_numeric)),
                           'type'=rep("Treatment window",length(treatment_days_numeric)))
  
  
  
  #we initialize the plot
  p <- ggplot()
  
  #The order of input matters a lot to ggplot and what happens behind the scenes is factors
  #get converted to numerics with appropriate labels, that you can then add continuous data to
  #you *cannot* start on a continuous scale and then add factors to it, so if we want to add a highlight
  #to the layer that's the furthest back we will need to initialize the plot as 
  #so we'll add an series of transparent dotplots with the appropriate x coordinates at an arbritary location on the y-axis
  p <- p + geom_point(data = unique(diversity_plotdata[,c("Day"),drop=F]),mapping = aes(x=Day,y=mean(diversity_plotdata[["Diversity"]])),alpha=0)
  
  #we start by highlighting the treatment days and add a new fill scale to include the boxplot fills
  p <- p + geom_ribbon(data = ribbon_dat,mapping = aes(x=treatment_days,ymax=ymax,ymin=ymin,fill=type),alpha=highlight_alpha,show.legend = F) + scale_fill_manual('Highlight',values = c('Treatment window'=treatment_window_highlight_color))#we add the background highlight first
  p <- p + ggnewscale::new_scale_fill()
  
  
  #we add borders around the window itself to clearly differentiate
  p <- p + geom_vline(xintercept = min(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color) + geom_vline(xintercept = max(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color)
  
  #first we add the boxplots but we don't add the outliers since we want to add those as colored points
  p <- p + geom_boxplot(diversity_plotdata,mapping = aes(x=Day,y=Diversity,fill=Group),width=boxplot_width,outlier.color = NA,linewidth = 0.3)
  
  #then we add the relevant points for the outliers and we do not add a legend because frankly why would we
  p <- p + geom_point(data = dotdat,mapping = aes(x=day_as_number,y=Diversity,color=Group,shape=Group),size=1.35,show.legend = F)
  p <- p + scale_shape_manual(values = cohort_pointshapes)
  
  #we then scale the color of the points and the fill of the boxplots to be consistent with our color-scheme
  p <- p  + scale_fill_manual('Dosage group',values = cohort_colors,breaks = names(cohort_colors),labels=cohort_scale_labels)+scale_color_manual('Dosage group',values = cohort_colors,breaks = names(cohort_colors),labels=cohort_scale_labels)
  
  #We set the y axis title and set additional breaks
  p <- p + scale_y_continuous(name = y_axis_text,breaks = seq(floor(plot_bottom),ceiling(plot_top),y_breaks_scale),limits = c(plot_bottom,plot_top),expand = c(0,0))
  
  #we add an x-axis title
  p <- p + xlab('Day')
  
  #and set our theme
  p <- p + custom_theme
  
  p <- p + annotate("label",x=mean(ribbon_dat[["treatment_days"]]),y=(plot_top-plot_bottom)*0.95+plot_bottom,label='Treatment window',color=treatment_highligt_text_color,hjust=0.5,size=2.8)
  
  p <- p + ggsignif::geom_signif(stat='identity',data = diversity_comparisons,mapping = aes(x=group1_as_number,
                                                                                            xend=group2_as_number,
                                                                                            annotation=as.character(signif(adjusted_p,2)),group=paste0(day,group2),
                                                                                            y=y,
                                                                                            yend=y),textsize = 1)
  return(list('figure'=p,'pvals'=diversity_comparisons))
}

#we create a function to create this type of recovery figure in a way that 
#dodges points that are closer to each other
recovery_figure_maker_pointranges_better_approach <- function(recovery_PFU_data,LLOQ,cohort_colors,cohort_scale_labels,cohort_pointshapes,plot_bottom,plot_top,type,dodge_extent,treatment_days){
  #we'll create a bunch of breaks and then we'll only add names to specific ones
  #these are the breaks we want but not have a label
  breaks_to_hide <- unlist(lapply(seq(-10,10),FUN = function(x){
    seq(2,9)*10^x
  }))
  #these are the breaks we wish to show
  breaks_to_show <- 10^seq(-10,10)
  
  #so we create empty break labels and proper break labels and combine them 
  labels <- c(to_tenth_power_labels(seq(-10,10)),rep('',length(breaks_to_hide)))
  breaks <- c(breaks_to_show,breaks_to_hide)
  
  #we make sure the names of the input data is standardized
  colnames(recovery_PFU_data) <- c('Subject.ID','Visit','Treatment','VALUE')
  
  #we will then convert the DAY X to just be a numeric X
  recovery_PFU_data[["Visit"]] <- as.numeric(gsub('DAY ','',recovery_PFU_data[["Visit"]]))
  
  #we standardize pre-screening data-points as -1
  recovery_PFU_data[["Visit"]][recovery_PFU_data[["Visit"]]%in%c(-1,-2)] <- -1
  
  #and we factorize the visit time since we do not especially care about the time-line in a numeric sense
  #as much as we care about the order
  recovery_PFU_data[["Visit"]] <- factor(recovery_PFU_data[["Visit"]])
  
  #We remove the SNIPR001 part of the dosage levels so we are just left with the dose level
  recovery_PFU_data[["Treatment"]] <- gsub(' SNIPR001','',recovery_PFU_data[["Treatment"]])
  
  #we'll turn the visit levels into numerics so we can treat the data as numeric and
  #manually dodge the point ranges and the path points later on
  recovery_PFU_data <- cbind(recovery_PFU_data,'Visit_NUMERIC'=as.numeric(recovery_PFU_data[["Visit"]]))
  
  #we create the breaks for the x-axis
  x_breaks <- seq(1,length(levels(recovery_PFU_data[["Visit"]])))
  #and then the labels that those breaks need to have
  x_labels <- levels(recovery_PFU_data[["Visit"]])
  x_labels[x_labels=='-1'] <- 'Baseline'
  
  #we'll add 1 to all the baseline values then log-transform the data at base 10
  recovery_PFU_data[["VALUE"]] <- log(recovery_PFU_data[["VALUE"]] + 1,base = 10)
  
  #we get the mean and sds of the data
  plotdat <- recovery_PFU_data %>% group_by(Visit_NUMERIC,Treatment) %>% summarise('mean'=mean(VALUE,na.rm=T),'sd'=sd(VALUE,na.rm=T))
  
  #and we'll change the order of the treatment groups
  plotdat[["Treatment"]] <- factor(plotdat[["Treatment"]],levels = c('Placebo','Cohort 1','Cohort 2','Cohort 3'))
  plotdat <- plotdat[order(plotdat[["Treatment"]]),]
  
  #we get the upper and lower bounds
  plotdat <- cbind(plotdat,'upper'=plotdat[["mean"]]+plotdat[["sd"]],'lower'=plotdat[["mean"]]-plotdat[["sd"]])
  
  #and then we offset as appropriate
  new_plotdat <- c()
  for(day in unique(plotdat[["Visit_NUMERIC"]])){
    offsets <- intelligently_dodge_one_timepoint(point_data = plotdat[plotdat[["Visit_NUMERIC"]]==day,c("Treatment","mean","upper","lower")],close_dodge = dodge_extent,whisker_point_dodge =dodge_extent*0.7,whisker_dodge = dodge_extent*0.27)
    new_plotdat <- rbind(new_plotdat,cbind(plotdat[plotdat[["Visit_NUMERIC"]]==day,],'offset'=offsets))
  }
  
  plotdat <- new_plotdat
  plotdat <- cbind(plotdat,'x_coordinate'=plotdat[["Visit_NUMERIC"]]+plotdat[["offset"]])
  
  #we'll get the numeric treatment days
  treatment_days_numeric <- seq(1,length(levels(recovery_PFU_data[["Visit"]])))[levels(recovery_PFU_data[["Visit"]])%in%treatment_days]
  
  #we'll get the lowest and highest x-coordinate in the treatment days and then add a bit to help 
  x_coord_range <- c(min(plotdat[plotdat[["Visit_NUMERIC"]]%in%treatment_days_numeric,"x_coordinate"])-0.15,
                     max(plotdat[plotdat[["Visit_NUMERIC"]]%in%treatment_days_numeric,"x_coordinate"])+0.15)
  
  #we'll create the ribbon data
  ribbon_dat <- data.frame('treatment_days'=x_coord_range,
                           'ymax'=rep(plot_top,length(x_coord_range)),
                           'ymin'=rep(plot_bottom,length(x_coord_range)),
                           'type'=rep("Treatment window",length(x_coord_range)))
  
  
  #first we initiate the figure
  p <- ggplot()
  #this figure is quite dense so we'll remove the gridlines that indicate logs while keeping the ticks
  #and just re-introduce some of them
  p <- p + geom_hline(yintercept = seq(-10,10),linewidth=0.15 ,color='lightgrey')+geom_vline(xintercept = x_breaks,linewidth=0.15,color='lightgrey')
  #we'll add the LLOQ line as well as a more well-defined zero line
  p <- p + geom_hline(yintercept = log(LLOQ,10),linetype='dashed') + geom_hline(yintercept = log(1,10),linetype='solid',alpha=0.4)
  #then we add the background highlight to identify the treatment window
  p <- p + geom_ribbon(data = ribbon_dat,mapping = aes(x=treatment_days,ymax=ymax,ymin=ymin,fill=type),alpha=highlight_alpha,show.legend=F) + scale_fill_manual('Highlight',values = c('Treatment window'=treatment_window_highlight_color))
  #we add borders around the window itself to clearly differentiate
  p <- p + geom_vline(xintercept = min(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color) + geom_vline(xintercept = max(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color)
  #we then add the lines to connect the different points with 
  p <- p + geom_path(data = plotdat,mapping = aes(x=x_coordinate,y=mean,color=Treatment,group=Treatment))
  #we'll add the points with the standard deviations to the same positions as the path points and add
  #shapes to help those who will print out in a non-color printer who like to print things
  p <- p + geom_pointrange(data = plotdat,mapping = aes(x=x_coordinate,y=mean,ymax=mean+sd,ymin=mean-sd,color=Treatment,group=Treatment,shape=Treatment),size=0.7)
  #we'll define the colors we use
  p <- p + scale_color_manual("Dosage",values = cohort_colors,breaks = names(cohort_colors),labels=cohort_scale_labels) + scale_x_continuous(breaks=x_breaks,labels=x_labels)
  #we'll set the shapes manually
  p <- p + scale_shape_manual('Dosage',values = cohort_pointshapes,breaks = names(cohort_colors),labels=cohort_scale_labels)
  #We'll set the y-axis labels, titles, and add a second scale to label the LLOQ on the right side of the graph instead of the left
  p <- p + scale_y_continuous(paste0("Mean ",type," SNIPR001 titer +/- SD (PFU/g)"),breaks=log(breaks,10),labels = labels,limits = limits,expand = c(0,0),sec.axis = sec_axis(trans = ~.,breaks = log(LLOQ,10),labels = paste0('LLOQ = ',LLOQ)))
  #and we set the x-axis title
  p <- p + xlab('Day')
  #we'll set the themes
  p <- p + custom_theme + theme(panel.grid.major = element_blank())
  #and return the figure
  
  return(p)
}

#we'll define the days where people are undergoing treatment
treatment_days <- c('1','3','5','7')

#and we'll define the labels to use for the legends 
cohort_scale_labels=c('Placebo'=expression('Placebo'),
                      'Cohort 1'=expression(10^8~'SNIPR001 PFU'),
                      'Cohort 2'=expression(10^10~'SNIPR001 PFU'),
                      'Cohort 3'=expression(10^12~'SNIPR001 PFU'))


#we define the relevant folders#
base_folder <- "your/folder/here/"
data_folder <- paste0(base_folder,"data/")
figures_folder <- paste0(base_folder,"figures/")
tables_folder <- paste0(base_folder,"tables/")


#We define our color scheme and other aesthetic coices
cohort_colors <- c('Placebo'="#EE9A00",'Cohort 1'="#004773",'Cohort 2'="#03858a",'Cohort 3'="#aae6c1")
treatment_window_highlight_color='#b3e2e8'
treatment_highligt_text_color="#45B3E6"
highlight_alpha=0.2               
cohort_pointshapes <- c('Placebo'=15,'Cohort 1'=16,'Cohort 2'=17,'Cohort 3'=18)



#We set a consistent ggplot theme in our figures called custom theme
text_size=8
custom_theme <- theme_bw()+#starts black white theme as starting point
  theme(axis.title.y=element_text(hjust=0.5,vjust=0.5,size=text_size),#Sets y axis title to be horizontal and appropriate size
        axis.title.x=element_text(size=text_size,angle=0),#sets x axis title
        axis.text.x = element_text(size=text_size),axis.text.y = element_text(size=text_size),#sets x axis text size
        legend.text = element_text(size=text_size),#sets legend text size
        panel.grid.minor.y=element_blank(),panel.grid.minor.x=element_blank(),#removes the minor panel grid for a cleaner look
        legend.title = element_text(size=text_size),#sets legend title
        legend.key.size = unit(0.5,'cm'),#sets the sie of the legend keys to be a half centimeter
        strip.background.x=element_rect(fill='#FFFFFF'),strip.background.y = element_rect(fill="#FFFFFF"),#sets facet backgrouns to white
        strip.text.x=element_text(size=text_size),strip.text.y = element_text(size=text_size),#sets facet title square backgrounds to be white
        panel.grid.major = element_line(linewidth = 0.1,color='#c9c9c9'),#sets the grid lines to be c pronounced
        axis.ticks = element_line(linewidth=0.1),#and sets the ticks to match the thickness
        legend.position='bottom')

multifigure_label_size <- 12

#we set the LLOQs for the SNIPR001 recovery assays as written in the protocols
stool_SNIPR001_recovery_LLOQ <- 4*10^2
plasma_SNIPR001_recovery_LLOQ <- 3*10^1
urine_SNIPR001_recovery_LLOQ <- 5*10^0



########################################
#Figure 1 drop-put stats and categories#
########################################
#first we'll get the reasons for screening failure
#we'll create the summary stuff
reasons_for_screen_fail <- openxlsx::read.xlsx(paste0(data_folder,'reasons_for_screen_failiure.xlsx'))
reasons_for_screen_fail <- reasons_for_screen_fail[,c("Subject.#","Status","Standardized.reason")]

#we'll get the categories we used for the paper
more_generalized_exclusion_reasions <-c(
  "Low E. coli" = "low E. coli",
  "Poor kidney function" ="not generally healthy",
  "Low bilirubin" ="not generally healthy",
  "Elevated BP" ="not generally healthy",
  "Low HR"="not generally healthy",
  "Current use of medication" ="used disallowed probiotics or medication",
  "High BMI"="not generally healthy",
  "High blood sugar"="not generally healthy",
  "Elevated total bilirubin"="not generally healthy",
  "Stool sample never recieved" ="sample issues",
  "Drug screen failure" ="used disallowed probiotics or medication",
  "Reseve"="reserve",
  "Elevated HR" ="not generally healthy",
  "Blood draw issues" ="sample issues",
  "No susceptible E. coli"="no susceptible E. coli",
  "Stool sample issues" ="sample issues",
  "High bristol stool scale"="not generally healthy",
  "Recent blood donation" ="unable to follow study protocol",
  "Allergy" ="not generally healthy",
  "Abnormal ECG"="not generally healthy",
  "High respiratory rate" ="not generally healthy",
  "Positive COVID test" ="other",
  "Unable to follow protocol" ="unable to follow study protocol",
  "Reserve" ="reserve",
  "Enrollment closed" ="enrolled outside trial window",
  "High amylase"="not generally healthy",
  "Proteinuria" ="not generally healthy",
  "Low HCT/HGB ratio" ="not generally healthy",
  "Blood in stool"="not generally healthy",
  "Pregnancy" ="other",
  "Urine sample issues" ="sample issues",
  "Elevated platelets"="not generally healthy",
  "Hepatitis" ="not generally healthy",
  "Recent participation in other drug trial"="other",
  "Nicotine use"="used disallowed probiotics or medication",
  "Recent probiotic use"="used disallowed probiotics or medication",
  "Withdrew consent"="withdrew consent",
  "ECG issues"="not generally healthy",
  "Significant medical findings"="not generally healthy",
  "History of relevant illness" ="not generally healthy",
  "Poor liver function" ="not generally healthy",
  "Not generally healthy" ="not generally healthy",
  "Heart related issues"="not generally healthy",
  "Poor pancreatic function"="not generally healthy",
  "Anemia"="not generally healthy",
  "High lipase"="not generally healthy",
  "Dosed"="Dosed",
  "Placebo"="Placebo"
)


#we'll add the more generalized reasons to the reasons for screen fail
reasons_for_screen_fail <- cbind(reasons_for_screen_fail,'Generalized_reason'=more_generalized_exclusion_reasions[reasons_for_screen_fail[,"Standardized.reason"]])

#since some people failed the screening for mulitple reaons we'll check if the any subjects had failed for multiple different generalized reasons
unique_subjects_reasons_ID_combos <- unique(reasons_for_screen_fail[,c("Subject.#","Status","Generalized_reason")])

#if anyone has been found to have failed/passed for multiple reasons we'll throw an error since it won't fit in the format of the table
if(max(table(unique_subjects_reasons_ID_combos[["Subject.#"]]))!=1){
  stop('Error! Some people failed for multiple reasons that are not in the same generalized reason')
}


#we'll write a sorted list of the reasons why people didn't pass or 
output_table <- table(unique_subjects_reasons_ID_combos[unique_subjects_reasons_ID_combos[["Status"]]!='Screen Pass',"Generalized_reason"])
output_table <- as.data.frame(output_table)
colnames(output_table) <- c('Categorized reason for failure','Number of failures')
write.xlsx(x = output_table,file = paste0(tables_folder,'Figure_1_categorized_reasons_for_failiure.xlsx'))



###################################################
#Table 1 summary demographics and plaquing results#
###################################################
#we load the relevant demograhics
relevant_demographics <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,'relevant_demographics.xlsx'))

#first we'll get the basic demographics data - male, female, age, and BMI out of the way
#first for the individual cohorts
demographics_data_distribution_cohort <- relevant_demographics %>% group_by(Actual.Treatment.for.Period) %>% summarise('n_male'=sum(Sex=="M"),
                                                                                                                          'n_female'=sum(Sex=="F"),
                                                                                                                          'mean_age_sd'=paste0(format(round(mean(Age),digits = 1)),'+-',format(round(sd(Age),digits = 1))),
                                                                                                                          'mean_BI_sd'=paste0(format(round(mean(`Body.Mass.Index.at.Baseline.(kg/m2)`),digits = 1)),'+-',format(round(sd(`Body.Mass.Index.at.Baseline.(kg/m2)`),digits = 1))))


#then the study as a whole
demographics_data_distribution_total <- relevant_demographics %>% summarise('n_male'=sum(Sex=="M"),
                                                                            'n_female'=sum(Sex=="F"),
                                                                            'mean_age_sd'=paste0(format(round(mean(Age),digits = 1)),'+-',format(round(sd(Age),digits = 1))),
                                                                            'mean_BI_sd'=paste0(format(round(mean(`Body.Mass.Index.at.Baseline.(kg/m2)`),digits = 1)),'+-',format(round(sd(`Body.Mass.Index.at.Baseline.(kg/m2)`),digits = 1))))


#and we gather ethnic information
ethnicity_distribution_cohort <- relevant_demographics %>% group_by(Actual.Treatment.for.Period,Ethnicity) %>% summarise('n'=length(Ethnicity)) %>% pivot_wider(id_cols = Actual.Treatment.for.Period,names_from = Ethnicity,values_from = n)
ethnicity_distribution_cohort[is.na(ethnicity_distribution_cohort)] <- 0


ethnicity_distribution_total <- cbind('All'='All',relevant_demographics) %>% group_by(All,Ethnicity) %>% summarise('n'=length(Ethnicity)) %>% pivot_wider(id_cols = "All",names_from = Ethnicity,values_from = n)
colnames(ethnicity_distribution_total)[1] <- colnames(ethnicity_distribution_cohort)[1]

ethnicity_distribution <- rbind(ethnicity_distribution_cohort,ethnicity_distribution_total)

#and race informmation
race_distribution_cohort <- relevant_demographics %>% group_by(Actual.Treatment.for.Period,Race) %>% summarise('n'=length(Race)) %>% pivot_wider(id_cols = Actual.Treatment.for.Period,names_from = Race,values_from = n)
race_distribution_cohort[is.na(race_distribution_cohort)] <- 0


race_distribution_total <- cbind('All'='All',relevant_demographics) %>% group_by(All,Race) %>% summarise('n'=length(Race)) %>% pivot_wider(id_cols = "All",names_from = Race,values_from = n)
colnames(race_distribution_total)[1] <- colnames(race_distribution_cohort)[1]

race_distribution <- rbind(race_distribution_cohort,race_distribution_total)

#we'll collapse the data together
age_sex_BMI_demographics_data <- rbind(demographics_data_distribution_cohort,unlist(c('all',unname(demographics_data_distribution_total))))


#we add the race and ethnicity info to the other demographics data
age_sex_BMI_demographics_data <- cbind(age_sex_BMI_demographics_data,
                                       race_distribution[,seq(2,ncol(ethnicity_distribution))],
                                       ethnicity_distribution[,seq(2,ncol(ethnicity_distribution))])


#we'll read the plaquing susceptibility data
all_quantitative_culture_data <- read.xlsx(xlsxFile =  paste0(data_folder,'quantiative ecoli and SNIPR001 data.xlsx'))

#we'll get all the culturing data from screening
screening_plaquing_data <- all_quantitative_culture_data[grepl('Number of',all_quantitative_culture_data[,"Parameter"]) & all_quantitative_culture_data[["Visit_name"]]=='Screening',]

#we'll get the plaquing data seperated into cohorts
plaquing_data_cohorts <- screening_plaquing_data %>% group_by(Subject.ID) %>% mutate('n_norm'=Value/sum(Value)) %>% group_by(Treatment,Parameter) %>% summarise('mean'=mean(n_norm*10),'sd'=sd(n_norm*10))

#and then the data for all together
plaquing_data_all <- screening_plaquing_data %>% group_by(Subject.ID) %>% mutate('n_norm'=Value/sum(Value)) %>% group_by(Parameter) %>% summarise('mean'=mean(n_norm*10),'sd'=sd(n_norm*10))
plaquing_data_all <- cbind('Treatment'='All',plaquing_data_all)


#we'll combine the tables
summary_plaquing_data <- rbind(plaquing_data_cohorts,plaquing_data_all)
summary_plaquing_data <- cbind(summary_plaquing_data,'table_text'=paste0(round(summary_plaquing_data[["mean"]],2),'±',round(summary_plaquing_data[["sd"]],2)))

#we reshape and output
summary_plaquing_data <- acast(data = summary_plaquing_data,formula = Parameter~Treatment,value.var = "table_text")
summary_plaquing_data <- summary_plaquing_data[c("Number of Plaque","Number of Lysis Zone","Number of Blank"),c("Cohort 1 SNIPR001","Cohort 2 SNIPR001","Cohort 3 SNIPR001","Placebo","All")]

#we get the distribution of stool titers at baseline
stool_titers_grouped <- all_quantitative_culture_data %>% filter(Parameter=='Escherichia coli' & Visit_name=='Screening') %>% group_by(Treatment) %>% summarise('mean'=mean(log(Value+1,10)),'sd'=sd(log(Value+1,10)))
stool_titers_all <- all_quantitative_culture_data %>% filter(Parameter=='Escherichia coli' & Visit_name=='Screening') %>% summarise('mean'=mean(log(Value+1,10)),'sd'=sd(log(Value+1,10)))
stool_titers_all <- c('Treatment'='All',stool_titers_all)

stool_titers_screening <- rbind(stool_titers_grouped,stool_titers_all)

write.xlsx(x = stool_titers_screening,file = paste0(tables_folder,'General_info_ecoli_titers_at_screening.xlsx'))

write.xlsx(x = age_sex_BMI_demographics_data,file = paste0(tables_folder,'Table_1_demographics_table_info_for.xlsx'))

#we save the summary plaquing data
write.xlsx(x = summary_plaquing_data,file = paste0(tables_folder,'Table_1_summary_plaquing_data.xlsx'))


########################################
#Table 2 solicited adverse events table#
########################################

#We load the adverse event data
adverse_events_data <- openxlsx::read.xlsx(paste0(data_folder,'adverse_events_data.xlsx'))

#We isolate the solicited adverse events and group the data by
#individual and type of adverse event
#we intentionally do not group by day as well, since some individuals experience
#persistent symptoms i.e. syptoms that persist over multiple days.
#These should be counted as a single persistent/recourring event
adverse_events_export <- adverse_events_data %>% filter(Solicited.AE=='Y') %>% group_by(Actual.treatment,Body.System.or.Organ.Class,Dictionary.Derived.Term) %>% summarise('n_events'=sum(Solicited.AE=='Y')) %>% arrange(Actual.treatment)

write.xlsx(x = adverse_events_export,file = paste0(tables_folder,'Table_2_adverse_events_summary_for.xlsx'))



########################################################
#Figure 3 - Recovery from blood, urine and stool figure#
########################################################

#we'll load all the recovery data
quant_ecoli_dat <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,'quantiative ecoli and SNIPR001 data.xlsx'))

#In order to ensure that each row is treated as a unique value we'll assign an ID to each measure. This will alleviate some of the issues around inconsistent sampling
quant_ecoli_dat <- cbind(quant_ecoli_dat,'measure_ID'=seq(1,nrow(quant_ecoli_dat)))

#er separate the relevant datasets 
stool_PFU_data <- quant_ecoli_dat[,c("Subject.ID","Parameter","Treatment","Value","Visit","measure_ID")] %>% filter(Parameter=='Stool Average titer') %>% tidyr::pivot_wider(names_from = Parameter,values_from=Value)
urine_PFU_data <- quant_ecoli_dat[,c("Subject.ID","Parameter","Treatment","Value","Visit","measure_ID")] %>% filter(Parameter=='Urine Average titer') %>% tidyr::pivot_wider(names_from = Parameter,values_from=Value)
plasma_PFU_data <- quant_ecoli_dat[,c("Subject.ID","Parameter","Treatment","Value","Visit","measure_ID")] %>% filter(Parameter=='Plasma Average titer') %>% tidyr::pivot_wider(names_from = Parameter,values_from=Value)

#We would like the different figures to shaare a common y axis limits
#in order to do that we'll first compute the highes tand lowest y-axis values and
#then use those as limits for all the figures and add a little space to each edge
temp <- quant_ecoli_dat %>% filter(grepl('Average titer',Parameter)) %>% mutate('log10_aval'=log(Value+1,10)) %>% group_by(Parameter,Visit_name,Treatment)  %>% summarise('top_point'=max(mean(log10_aval)+sd(log10_aval)),'bottom_point'=min(mean(log10_aval)-sd(log10_aval)))
limits <- c(min(temp[["bottom_point"]])-0.5,max(temp[["top_point"]])+1.2)

#we'll set the extent to which we dodge the points
dodge_extent <- 0.20

#we'll create each of the figures for stool, plasma, and urine, respectively
stool_figure <- recovery_figure_maker_pointranges_better_approach(recovery_PFU_data  = stool_PFU_data[,c("Subject.ID","Visit","Treatment","Stool Average titer")],LLOQ = stool_SNIPR001_recovery_LLOQ,cohort_colors = cohort_colors,cohort_pointshapes = cohort_pointshapes,cohort_scale_labels = cohort_scale_labels,plot_bottom = limits[1],plot_top = limits[2],type = 'stool',dodge_extent=dodge_extent*0.7,treatment_days = treatment_days)
plasma_figure <- recovery_figure_maker_pointranges_better_approach(recovery_PFU_data = plasma_PFU_data[,c("Subject.ID","Visit","Treatment","Plasma Average titer")],LLOQ =plasma_SNIPR001_recovery_LLOQ,cohort_colors = cohort_colors,cohort_pointshapes = cohort_pointshapes,cohort_scale_labels = cohort_scale_labels,plot_bottom = limits[1],plot_top = limits[2],type = 'plasma',dodge_extent=dodge_extent,treatment_days = treatment_days)
urine_figure <- recovery_figure_maker_pointranges_better_approach(recovery_PFU_data = urine_PFU_data[,c("Subject.ID","Visit","Treatment","Urine Average titer")],LLOQ =urine_SNIPR001_recovery_LLOQ,cohort_colors = cohort_colors,cohort_pointshapes = cohort_pointshapes,cohort_scale_labels = cohort_scale_labels,plot_bottom = limits[1],plot_top = limits[2],type = 'urine',dodge_extent=dodge_extent,treatment_days = treatment_days)

#we'll need to get the coordinates for the treatment window
treatment_window_data <- stool_figure[["layers"]][[5]][["data"]]
stool_figure <- stool_figure + annotate("label",x=mean(treatment_window_data[["treatment_days"]]),y=diff(limits)*0.95+min(limits),label='Treatment window',color=treatment_highligt_text_color,size=2.6)

#we'll manipulate the top figure even further introducing the legend into the figure itself
stool_figure <- stool_figure + theme(legend.position = 'none')

#these will be placed custom. No way I'm figuring out a way to solve the general case for this one lmao
stool_figure <- stool_figure + annotate(geom = 'label',x = 4.5,y = 4,label=expression(10^8~' PFU'),hjust=0.5,vjust=0,color=cohort_colors['Cohort 1'],size=2.6)
stool_figure <- stool_figure + annotate(geom = 'label',x = 4.5,y = 6.5,label=expression(10^10~' PFU'),hjust=0.5,vjust=0,color=cohort_colors['Cohort 2'],size=2.6)
stool_figure <- stool_figure + annotate(geom = 'label',x = 4.5,y = 9.3,label=expression(10^12~' PFU'),hjust=0.5,vjust=0,color=cohort_colors['Cohort 3'],size=2.6)
stool_figure <- stool_figure + annotate(geom = 'label',x = 4.5,y = 0.4,label="Placebo",hjust=0.5,vjust=0,color=cohort_colors['Placebo'],size=2.6)

#we want a figure that uses the stool recovery figure as a main figure, and then two sub-figures below. These all share the same
#legend so we'll just include the legend for the stool figure

#we'll assemble the bottom plots first, which is plasma and urine recovery without legends
bottom_plots <- cowplot::plot_grid(plotlist = list(plasma_figure+theme(legend.position = 'none'),urine_figure+theme(legend.position = 'none')),align = 'hv',axis = 'trbl',labels = c('B','C'),label_size = multifigure_label_size)

#we'll add the main stool figure
top_plot <- cowplot::plot_grid(plotlist = list(stool_figure),labels = c('A'),label_size = multifigure_label_size)

#then we merge the two levels in a single figure
figure <- cowplot::plot_grid(plotlist = list(top_plot,bottom_plots),align = 'hv',axis = 'trbl',ncol = 1)

#then we save the figure of all the recovery figure
ggsave(plot = figure,filename = paste0(figures_folder,'Figure_3_recovery_illustration.pdf'),width = 8,height = 6)


#########################################
#Figure 2A E. coli recovery by LS means#
#########################################
#we load Medpace's LS means data
lsmeans_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,'lmm_predictions_full_dataset.xlsx'))

#we'll factorize the days since we are not particularly interested in the numeric distance betwen
#each measurement as much as the order
lsmeans_data[["Day"]] <- factor(lsmeans_data[["Day"]])

#we set the position dodge width
dodge_width=0.35

#we get the ranges for the treatment highlight
plot_top <- max(lsmeans_data[["upper.95.percent.CI"]])+0.1
plot_bottom <- min(lsmeans_data[["lower.95.percent.CI"]])-0.1

#For some formatting stuff we'll set the x-axis data to be numeric and then
#set the x-axis labels and breaks to be what we use before numeric conversion
x_labels <- levels(lsmeans_data[["Day"]])
x_breaks <- seq(1,length(x_labels))
lsmeans_data[["Day_numeric"]] <- as.numeric(lsmeans_data[["Day"]])

#we get the numeric values of the treatment days
numeric_treament_days <- x_breaks[x_labels%in%treatment_days]

#we sort the cohorts
lsmeans_data[["Cohort"]] <- factor(lsmeans_data[["Cohort"]],levels = c('Placebo','Cohort 1','Cohort 2','Cohort 3'))

#then we extend the x-axis treatment highlight width a bit to account for the extent of the position_dodge we'll put in later
numeric_treament_days[1] <- numeric_treament_days[1]-(dodge_width/2)*1.2
numeric_treament_days[length(numeric_treament_days)] <- numeric_treament_days[length(numeric_treament_days)]+(dodge_width/2)*1.2

#we compile the ribbon data
ribbon_dat = data.frame('treatment_days'=numeric_treament_days,
                        'ymax'=rep(plot_top,length(treatment_days)),
                        'ymin'=rep(plot_bottom,length(treatment_days)),
                        'type'='Treatment window')

#we'll create a bunch of breaks and then we'll only add names to specific ones
#these are the breaks we want but not have a label
breaks_to_hide <- unlist(lapply(seq(-10,10),FUN = function(x){
  seq(2,9)*10^x
}))
#these are the breaks we wish to show
breaks_to_show <- 10^seq(-10,10)

#so we create empty break labels and proper break labels and combine them 
labels <- c(to_tenth_power_labels(seq(-10,10)),rep('',length(breaks_to_hide)))
breaks <- c(breaks_to_show,breaks_to_hide)

#we initialize the figure
p <- ggplot()

# In order to highlight the no-change condition we'll highlight the x-axis of no change 10^0
p <- p + geom_hline(yintercept = 0,linetype='dashed') #and we'll add a line to indicate zero change

#we start by highlighting the treatment days
p <- p + geom_ribbon(data = ribbon_dat,mapping = aes(x=treatment_days,ymax=ymax,ymin=ymin,fill=type),alpha=highlight_alpha,show.legend = F) + scale_fill_manual('Highlight',values = c('Treatment window'=treatment_window_highlight_color))#we add the background highlight first

#we add borders around the window itself to clearly differentiate
p <- p + geom_vline(xintercept = min(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color) + geom_vline(xintercept = max(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color)

#we'll add each LS means point with their 95% CI using point ranges as 95% CIs
p <- p + geom_pointrange(data=lsmeans_data,mapping = aes(x=Day_numeric,color=Cohort,y=lsmeans,group=Cohort,ymax=upper.95.percent.CI ,ymin=lower.95.percent.CI,shape=Cohort),linewidth=0.7,size=0.6,position = position_dodge2(width = dodge_width)) #then we add the points with 95% CI

#we'll set the shapes to be consistent with previous shapes 
p <- p + scale_shape_manual('Dosage',values = cohort_pointshapes,labels=cohort_scale_labels)

#we'll set the colors to be consisten with previous colors
p <- p + scale_color_manual("Dosage",values = cohort_colors,labels=cohort_scale_labels) #then we change the colors and set the theme

#we se the x and y axis titles
p <- p + xlab('Day')+ylab(substitute(paste('LS means prediction +/- 95% CI ',italic('E. coli'), ' count change from day -1 (CFU/g)',sep = ' '))) #we set the axis labels

#we'll scale the y axis to include log labels as the regression was carried out in log-scale and limit the y-axis to the range highlighted by the treatment highlight
p <- p + scale_y_continuous(breaks = log(breaks,10),labels = labels,limits = c(plot_bottom,plot_top),expand = c(0,0))#we set the y-axis axis labels and set the y-lim to be exactly between the ribbon start and end

#we'll scale the x-axis to include the initial labels before numeric conversion
p <- p + scale_x_continuous(breaks = x_breaks,labels = x_labels) # we also set the x axis breaks and re-introduce the original break labels

#we'll add a label that indicates the treatment window here 
p <- p + annotate("label",x=mean(ribbon_dat[["treatment_days"]]),y=(plot_top-plot_bottom)*0.95+plot_bottom,label='Treatment window',color=treatment_highligt_text_color,hjust=0.5,size=2.8)

#we'll set the theme
p_lm_figure_187_days <- p + custom_theme


#we'll compute the differences in means
difference_data <- lsmeans_data %>% group_by(Day) %>% summarise('Cohort'=Cohort,
                                                                'mean_diff'=lsmeans-lsmeans[Cohort=='Placebo'],
                                                                'SE_diff'=sqrt(SE^2+SE[Cohort=='Placebo']^2)) %>% 
                                                      filter(Cohort!='Placebo')



#next we'll calculate the chance of chis number being different from zero
difference_data[['two_tailed_pval_not_corrected']] <- 2*pnorm(-abs(difference_data[["mean_diff"]]/difference_data[["SE_diff"]]))

#we'll FDR correct
difference_data <- difference_data %>% group_by(Cohort) %>% mutate('FDR_corrected_Pval'=p.adjust(two_tailed_pval_not_corrected,method='fdr'))

colnames(difference_data) <- c('Day','Cohort compared to placebo','Difference of means','standard error of difference','Two-tailed P-value','FDR corrected P-value')

#we'll output this data
openxlsx::write.xlsx(x = difference_data,file = paste0(tables_folder,'Supplementary_table_2A_lmer_significance_data_full_dataset.xlsx'))

#################################
#Figure 2B - Beta diversity plot#
#################################
#we read the beta diversities comparing each sample to each other sample using unique sample IDs
beta_divs <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,'braycurtis_matrix.all.xlsx'),rowNames = T)

#we read the metadata of each sample
sample_metadata <- openxlsx::read.xlsx(xlsxFile =paste0(data_folder,'sample_metadata.xlsx'))


#we'll reformat a format that's easier for us to interact with
beta_div_melt <- reshape2::melt(as.matrix(beta_divs))
colnames(beta_div_melt) <- c('sample1_id','sample2_id','beta_div')

#To make interacting a bit easier we set the rownames of the meta data to be the sample ID
rownames(sample_metadata) <- sample_metadata[["unblinded_sample_id"]]

#we get the relevant metadata for the first sample in the beta diversity
sample1_append <- sample_metadata[beta_div_melt[["sample1_id"]],c("unblinded_subject_id","Day","cohort","treatment_group")]
colnames(sample1_append) <- paste0('sample1_',colnames(sample1_append))

#we get the relevant metadata for the second sample in the beta diversity
sample2_append <- sample_metadata[beta_div_melt[["sample2_id"]],c("unblinded_subject_id","Day","cohort","treatment_group")]
colnames(sample2_append) <- paste0('sample2_',colnames(sample2_append))

#then we append the relevant metadata to the beta diversities
beta_div_melt <- cbind(beta_div_melt,sample1_append,sample2_append)

#first we'll make sure we only compare a person to their previous time-point and that previous time-point is -1
beta_div_melt <- beta_div_melt[beta_div_melt[["sample1_unblinded_subject_id"]]==beta_div_melt[["sample2_unblinded_subject_id"]] & #compare samples to samples from own person
                                 beta_div_melt[["sample1_Day"]]=='-1' & #and only those where the first sample is the screening sample
                                 beta_div_melt[["sample2_Day"]]>0, #and we don't compare to samples before the screening day
] 

#we'll set the day as a factor to help the boxplot
beta_div_melt[["sample2_Day"]] <- factor(beta_div_melt[["sample2_Day"]],levels = sort(unique(beta_div_melt[["sample2_Day"]]),decreasing = F))

#We'll change the cohort to reflect the treatment recieved, not the timing of the dosign
beta_div_melt[beta_div_melt[["sample2_treatment_group"]]=='Placebo',c("sample2_cohort","sample1_cohort")] <- 'Placebo'

#And then we capitalize the cohort names
beta_div_melt[["sample1_cohort"]] <- gsub('cohort','Cohort ',beta_div_melt[["sample1_cohort"]])

#we'll order things properly
beta_div_melt[["sample1_cohort"]] <- factor(beta_div_melt[["sample1_cohort"]],levels = names(cohort_colors))


#we format and generate the plot
diversity_plotdata <- beta_div_melt[,c("sample2_Day","sample1_cohort","beta_div")]
colnames(diversity_plotdata) <- c('Day','Group','Diversity')

diversity_returns <- diversity_figure_generator(diversity_plotdata = diversity_plotdata,
                           y_axis_text = "Bray-Curtis dissimilarity to pre-dosage time-point",
                           treatment_days = treatment_days,
                           cohort_pointshapes = cohort_pointshapes,
                           cohort_colors = cohort_colors,
                           custom_theme = custom_theme,
                           cohort_scale_labels = cohort_scale_labels,
                           treatment_window_highlight_color = treatment_window_highlight_color,
                           highlight_alpha = highlight_alpha,y_breaks_scale = 0.1,comparison_type = 'two.sided')

p_beta_diversity <- diversity_returns[["figure"]]

##########
#Figure 2#
##########
#we'll make a unifying figure legend which would just be the colors Easiest way
#to do that is to just make a heatmap and steal the legend from that
legend_plotdat <- data.frame('type'=names(cohort_colors))
legend_figure <- ggplot(data = legend_plotdat,mapping = aes(x=type,y=type,fill=type))+geom_tile()+scale_fill_manual('Dosage',values = cohort_colors,breaks = names(cohort_colors),labels=cohort_scale_labels)+theme(legend.position = 'bottom')
unifying_legend <- get_legend(legend_figure)

#we'll combine the two figures
plotlist <- list(p_beta_diversity+theme(legend.position = 'none'),
                 p_lm_figure_187_days+theme(legend.position = 'none'),
                 unifying_legend)

comboplot <- cowplot::plot_grid(plotlist = plotlist,ncol = 1,rel_heights = c(1.1,1,0.1),align='hv',axis = 'trbl')

ggsave(plot = comboplot,filename = paste0(figures_folder,'Figure_2_e_coli_reduction_and_beta_div.pdf'),height = 7.5,width = 8)

################################################
#Supplementary figure 3 - Alpha diversity plots#
################################################
#and we'll do the alpha diversity change figure
alpha_diversities <- read.xlsx(xlsxFile = paste0(data_folder,'alpha_diversities.xlsx'))

#we'll only include the IDs who were measured at every the pre-dose timepoint
screening_subjects <- alpha_diversities[["unblinded_subject_id"]][alpha_diversities[["Day"]]==-1]
alpha_diversities <- alpha_diversities[alpha_diversities[["unblinded_subject_id"]]%in%screening_subjects,]


#we'll compute the alpha diversities relative to pre-treatment for each individual
alpha_diversity_plotdat <- alpha_diversities %>% group_by(unblinded_subject_id) %>% reframe('Day'=Day,
                                                                                            'Treatment'=Treatment,
                                                                                            'relative_shannon'=shannon_alpha-shannon_alpha[Day==-1],
                                                                                            'relative_Chao1'=chao1_alpha-chao1_alpha[Day==-1])

#we remove day -1 since it'll be 0s all around and -30 because it's not relevant
#to to how dosing changes alpha diversity
alpha_diversity_plotdat <- alpha_diversity_plotdat[!alpha_diversity_plotdat[["Day"]]%in%c('-1',"-30"),]

#we factorize the days
alpha_diversity_plotdat[["Day"]] <- factor(alpha_diversity_plotdat[["Day"]],levels = sort(unique(alpha_diversity_plotdat[["Day"]]),decreasing = F))

#we'll order things properly
alpha_diversity_plotdat[["Treatment"]] <- factor(alpha_diversity_plotdat[["Treatment"]],levels = names(cohort_colors))


diversity_plotdata <- alpha_diversity_plotdat[,c("Day","Treatment","relative_shannon")]
colnames(diversity_plotdata) <- c('Day','Group','Diversity')

p_shannon_diversity <- diversity_figure_generator(diversity_plotdata = diversity_plotdata,
                                               y_axis_text = "Relative Shannon diversity relative to pre-dose time-point",
                                               treatment_days = treatment_days,
                                               cohort_pointshapes = cohort_pointshapes,
                                               cohort_colors = cohort_colors,
                                               custom_theme = custom_theme,
                                               cohort_scale_labels = cohort_scale_labels,
                                               treatment_window_highlight_color = treatment_window_highlight_color,
                                               highlight_alpha = highlight_alpha,
                                               y_breaks_scale = 0.5,
                                               comparison_type='two.sided')

p_shannon_diversity <- p_shannon_diversity[["figure"]]

#next up we'll do the Chao1 index
diversity_plotdata <- alpha_diversity_plotdat[,c("Day","Treatment","relative_Chao1")]
colnames(diversity_plotdata) <- c('Day','Group','Diversity')


p_chao1_diversity <- diversity_figure_generator(diversity_plotdata = diversity_plotdata,
                                                y_axis_text = "Relative Chao1 diversity relative to pre-dose time-point",
                                                treatment_days = treatment_days,
                                                cohort_pointshapes = cohort_pointshapes,
                                                cohort_colors = cohort_colors,
                                                custom_theme = custom_theme,
                                                cohort_scale_labels = cohort_scale_labels,
                                                treatment_window_highlight_color = treatment_window_highlight_color,
                                                highlight_alpha = highlight_alpha,
                                                y_breaks_scale = 10,
                                                comparison_type='two.sided')

#due to the high scale of chao1 indeces this kind of breaks the y scale settings
#that we used before and makes it look kinda ugly. We fix it
p_chao1_diversity <- p_chao1_diversity[["figure"]]
p_chao1_diversity <- p_chao1_diversity + scale_y_continuous(name = "Relative Chao1 diversity relative to screening time-point",breaks = seq(-1000,1000,20),expand = c(0,0))

#we'll merge the legends for the combined plot so we first extract the legend
#and then remove it from both
common_legend <- get_legend(p_chao1_diversity)

#we'll put the two together
plotlist <- list(p_shannon_diversity+theme(legend.position = 'none'),
                 p_chao1_diversity+theme(legend.position = 'none'),
                 common_legend)


comboplot <- cowplot::plot_grid(plotlist = plotlist,align = 'hv',axis = 'trbl',ncol=1,rel_heights = c(1,1,.1))

#we save the figure
ggsave(plot = comboplot,filename = paste0(figures_folder,'Supplementary_figure_3_alpha_div.pdf'),width = 10,height = 8)


##############################################################################################
#Suppelementary figure 2 - Mixed linear model outcomes - the model based on the first 14 days#
##############################################################################################
#we load Medpace's LS means data
lsmeans_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,'lmm_predictions_up_to_day_14.xlsx'))

#we'll factorize the days since we are not particularly interested in the numeric distance betwen
#each measurement as much as the order
lsmeans_data[["Day"]] <- factor(lsmeans_data[["Day"]])

#we set the position dodge width
dodge_width=0.5

#we get the ranges for the treatment highlight
plot_top <- max(lsmeans_data[["upper.95.percent.CI"]])+0.1
plot_bottom <- min(lsmeans_data[["lower.95.percent.CI"]])-0.1

#For some formatting stuff we'll set the x-axis data to be numeric and then
#set the x-axis labels and breaks to be what we use before numeric conversion
x_labels <- levels(lsmeans_data[["Day"]])
x_breaks <- seq(1,length(x_labels))
lsmeans_data[["Day_numeric"]] <- as.numeric(lsmeans_data[["Day"]])

#we get the numeric values of the treatment days
numeric_treament_days <- x_breaks[x_labels%in%treatment_days]

#we sort the cohorts
lsmeans_data[["Cohort"]] <- factor(lsmeans_data[["Cohort"]],levels = c('Placebo','Cohort 1','Cohort 2','Cohort 3'))
lsmeans_data <- lsmeans_data[order(lsmeans_data[["Cohort"]]),]

#then we extend the x-axis treatment highlight width a bit to account for the extent of the position_dodge we'll put in later
numeric_treament_days[1] <- numeric_treament_days[1]-(dodge_width/2)*1.2
numeric_treament_days[length(numeric_treament_days)] <- numeric_treament_days[length(numeric_treament_days)]+(dodge_width/2)*1.2

#we compile the ribbon data
ribbon_dat = data.frame('treatment_days'=numeric_treament_days,
                        'ymax'=rep(plot_top,length(treatment_days)),
                        'ymin'=rep(plot_bottom,length(treatment_days)),
                        'type'='Treatment window')

#we'll create a bunch of breaks and then we'll only add names to specific ones
#these are the breaks we want but not have a label
breaks_to_hide <- unlist(lapply(seq(-10,10),FUN = function(x){
  seq(2,9)*10^x
}))
#these are the breaks we wish to show
breaks_to_show <- 10^seq(-10,10)

#so we create empty break labels and proper break labels and combine them 
labels <- c(to_tenth_power_labels(seq(-10,10)),rep('',length(breaks_to_hide)))
breaks <- c(breaks_to_show,breaks_to_hide)


#we initialize the figure
p <- ggplot()

# In order to highlight the no-change condition we'll highlight the x-axis of no change 10^0
p <- p + geom_hline(yintercept = 0,linetype='dashed') #and we'll add a line to indicate zero change

#we start by highlighting the treatment days
p <- p + geom_ribbon(data = ribbon_dat,mapping = aes(x=treatment_days,ymax=ymax,ymin=ymin,fill=type),alpha=highlight_alpha,show.legend = F) + scale_fill_manual('Highlight',values = c('Treatment window'=treatment_window_highlight_color))#we add the background highlight first

#we add borders around the window itself to clearly differentiate
p <- p + geom_vline(xintercept = min(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color) + geom_vline(xintercept = max(ribbon_dat[["treatment_days"]]),color=treatment_highligt_text_color)

#we'll add each LS means point with their 95% CI using point ranges as 95% CIs
p <- p + geom_pointrange(data=lsmeans_data,
                         mapping = aes(x=Day_numeric,color=Cohort,y=lsmeans,group=Cohort,ymax=upper.95.percent.CI ,ymin=lower.95.percent.CI ,shape=Cohort),
                         size=0.6,position = position_dodge2(width = dodge_width)) #then we add the points with 95% CI

#we'll set the shapes to be consistent with previous shapes 
p <- p + scale_shape_manual('Dosage',values = cohort_pointshapes,labels=cohort_scale_labels)

#we'll set the colors to be consisten with previous colors
p <- p + scale_color_manual("Dosage",values = cohort_colors,labels=cohort_scale_labels) #then we change the colors and set the theme

#we se the x and y axis titles
p <- p + xlab('Day')+ylab(substitute(paste('LS means prediction +/- 95% CI ',italic('E. coli'), ' \n counts change from day -1 (CFU/g)',sep = ' '))) #we set the axis labels

#we'll scale the y axis to include log labels as the regression was carried out in log-scale and limit the y-axis to the range highlighted by the treatment highlight
p <- p + scale_y_continuous(breaks = log(breaks,10),labels = labels,limits = c(plot_bottom,plot_top),expand = c(0,0))#we set the y-axis axis labels and set the y-lim to be exactly between the ribbon start and end

#we'll scale the x-axis to include the initial labels before numeric conversion
p <- p + scale_x_continuous(breaks = x_breaks,labels = x_labels) # we also set the x axis breaks and re-introduce the original break labels

#we'll set the theme
p <- p + custom_theme

p <- p + annotate("label",x=mean(ribbon_dat[["treatment_days"]]),y=(plot_top-plot_bottom)*0.95+plot_bottom,label='Treatment window',color=treatment_highligt_text_color,hjust=0.5,size=2.8)

#we save the figure
ggsave(plot = p,filename = paste0(figures_folder,'Supplementary_figure_2_lmer_up_to_day_14_model.pdf'),width = 8,height = 5)



#we'll compute the differences in means
difference_data <- lsmeans_data %>% group_by(Day) %>% summarise('Cohort'=Cohort,
                                                                'mean_diff'=lsmeans-lsmeans[Cohort=='Placebo'],
                                                                'SE_diff'=sqrt(SE^2+SE[Cohort=='Placebo']^2)) %>% 
  filter(Cohort!='Placebo')



#next we'll calculate the chance of chis number being different from zero
difference_data[['two_tailed_pval_not_corrected']] <- 2*pnorm(-abs(difference_data[["mean_diff"]]/difference_data[["SE_diff"]]))

#we'll FDR correct
difference_data <- difference_data %>% group_by(Cohort) %>% mutate('FDR_corrected_Pval'=p.adjust(two_tailed_pval_not_corrected,method='fdr'))

colnames(difference_data) <- c('Day','Cohort compared to placebo','Difference of means','standard error of difference','Two-tailed P-value','FDR corrected P-value')

#we'll output this data
openxlsx::write.xlsx(x = difference_data,file = paste0(tables_folder,'Supplementary_table_2B_lmer_significance_data_day_14_subset.xlsx'))


###################################
#Table 3 - Table of adverse events#
###################################
#We load the adverse event data up to Day 35
adverse_events_data <- openxlsx::read.xlsx(paste0(data_folder,'adverse_events_data.xlsx'))

#we select relevant columns and get the number of each
all_adverse_events <- adverse_events_data %>%  group_by(Actual.treatment,Solicited.AE,Analysis.Toxicity.Grade,Causality,Dictionary.Derived.Term,Outcome.of.Adverse.Event) %>% reframe('n'=length(Dictionary.Derived.Term)) %>% arrange(Actual.treatment,Solicited.AE,Causality,Analysis.Toxicity.Grade,Dictionary.Derived.Term)

#we reformat things slightly to spell things out for the article
all_adverse_events[["Actual.treatment"]] <- c("Cohort 1 SNIPR001"='10^8 BID',"Cohort 2 SNIPR001"='10^10 BID',"Cohort 3 SNIPR001"='10^12 BID','Placebo'='Placebo')[all_adverse_events[["Actual.treatment"]]]
all_adverse_events[["Causality"]] <- c('N'='No','Y'='Yes')[all_adverse_events[["Causality"]]]
all_adverse_events[["Solicited.AE"]] <- c('N'='No','Y'='Yes')[all_adverse_events[["Solicited.AE"]]]
colnames(all_adverse_events) <- gsub('[.]',' ',colnames(all_adverse_events))

#and we save it
openxlsx::write.xlsx(x = all_adverse_events,file = paste0(tables_folder,'Table_3_adverse_events_table_up_to_day_35.xlsx'))

#############################################
#Supplementary table 5 - AST/ALT ratio tests#
#############################################

bloodwork <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,'bloodwork data.xlsx'))


relevant_bloodwork <- bloodwork[bloodwork[["Parameter"]]%in%c("Aspartate Aminotransferase (U/L)","Alanine Aminotransferase (U/L)"),]


relevant_bloodwork <- relevant_bloodwork %>%  pivot_wider(names_from = Parameter,values_from = Value)


relevant_bloodwork[['ast_alt_ratio']] <- relevant_bloodwork[["Aspartate Aminotransferase (U/L)"]]/relevant_bloodwork[["Alanine Aminotransferase (U/L)"]]

all_groups <- unique(relevant_bloodwork[["Treatment"]])

placebo_group <- 'Placebo'
non_placebo_groups <- all_groups[!all_groups==placebo_group]

days <- unique(relevant_bloodwork[["Visit"]])


ast_alt_comparison_data <- data.frame()
for(day in days){
  placebo_data <- relevant_bloodwork[["ast_alt_ratio"]][relevant_bloodwork[["Treatment"]]==placebo_group & relevant_bloodwork[["Visit"]]==day]
  
  for(comparison in non_placebo_groups){
    group_data <- relevant_bloodwork[["ast_alt_ratio"]][relevant_bloodwork[["Treatment"]]==comparison & relevant_bloodwork[["Visit"]]==day]
    
    test_result <- wilcox.test(group_data,placebo_data)
    
    diff <- effsize::cliff.delta(group_data,placebo_data)
    
    ast_alt_comparison_data <- rbind(ast_alt_comparison_data,data.frame('group1'=placebo_group,'group2'=comparison,'day'=day,'pval'=test_result[["p.value"]],'cliff_delta'=diff[["estimate"]]))
  }
}


ast_alt_comparison_data[["pval"]] <- round(ast_alt_comparison_data[["pval"]],2)
ast_alt_comparison_data[["cliff_delta"]] <- round(ast_alt_comparison_data[["cliff_delta"]],2)

openxlsx::write.xlsx(x = ast_alt_comparison_data,file = paste0(tables_folder,'Supplementary_table_5_ast_alt_ratio_pval_day.xlsx'))




##########################################
#Supplementary table 4 - test of MIC data#
##########################################
#first we read the data
amr_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_folder,"amr_data.xlsx"))

#let's try post and pre-treatment
amr_data <- cbind(amr_data,'treated_status'=c('Untreated','Treated')[as.numeric(amr_data[["Day"]]>0)+1])

#we'll set the placebo guys to be untreated at all timepoints
amr_data[["treated_status"]][amr_data[["unblinded_subject_id"]]%in%unique(quant_ecoli_dat[["Subject.ID"]][quant_ecoli_dat[["Treatment"]]=='Placebo'])] <- 'Untreated'

amr_data[["treated_status"]] <- factor(amr_data[["treated_status"]],levels = c('Treated',"Untreated"))


#We'll see if the MIC has anything to do with the treatment status
amr_linear_model <- lme4::lmer(data = amr_data,formula = "mic~treated_status*Drug+(1|strain_id)")

emm <- emmeans::emmeans(amr_linear_model,~treated_status*Drug)
contrasts <- emmeans::contrast(emm, method = "pairwise", adjust = "none") %>% as.data.frame

#we extract the relevant contrasts
relevant_rows <- lapply(contrasts[["contrast"]],FUN = function(x){
  split_dat <- strsplit(x,split = ' [-] ')[[1]]
  if(sum(grepl('Treated',split_dat))==1){
    first_space1 <- gregexpr(pattern = '[ ]',text = split_dat[1])[[1]][1]
    first_space2 <- gregexpr(pattern = '[ ]',text = split_dat[2])[[1]][1]
    
    return(substr(split_dat[1],start = first_space1+1,nchar(split_dat[1]))==substr(split_dat[2],start = first_space2+1,nchar(split_dat[2])))
  }else{
    return(F)
  }
}) %>% unlist

relevant_contrasts <- contrasts[relevant_rows,]

#we'll get the one-way test
one_way_pval <- apply(relevant_contrasts,1,FUN = function(x){
  pt(q = (as.numeric(x["estimate"])/as.numeric(x["SE"])),df = as.numeric(x["df"]),lower.tail = FALSE)
})

relevant_contrasts[['one_way_pval_greater']]=one_way_pval

relevant_contrasts[['FDR_adjusted_pval']] <- p.adjust(relevant_contrasts[['one_way_pval_greater']],method = 'fdr')


#we'll export
openxlsx::write.xlsx(x = relevant_contrasts[,c("contrast","estimate","SE","one_way_pval_greater","FDR_adjusted_pval")],file = paste0(tables_folder,'Supplementary_table_4_abx_susceptibility_regression_for_supplementary.xlsx'))



