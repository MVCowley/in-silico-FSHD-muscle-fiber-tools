# This file is part of a set of shiny applications for the analysis of in silico
# FSHD muscle fibres.
# Copyright © 2023 C. R. S. Banerji
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

library(shiny)
library("tidyverse")
library('deSolve')
library(lattice)
library('ggplot2')
load("files_for_CA_MT_size_shiny_not_refined.rd")
server<-function(input,output,session){
  ##############
  #########define reative values which update as the automaton evolves
  rvs<-reactiveValues(iter=array("S",dim=c(mt_l,mt_c)),bar=setNames(c(mt_l*mt_c,0,0,0,0),c("S","E","I","R","D")),ret3=ret3,
                      counter=1,reset=0)
  #########output dynamically updating plot of automaton at the top
  output$plot1<-renderPlot({
    plot(plot_realisation(rvs$iter))})
  #########output dynamically updating histogram of the automaton underneath
  output$plot2<-renderPlot({
    ####organise
    temp4<-data.frame(col_anno[,1],rep(0,5),col_anno[,2])
    colnames(temp4)<-c("State","Values","Col")
    temp4$State<-factor(temp4$State,levels=col_anno[,1])
    which(is.na(match(col_anno[,1],names(rvs$bar))))->m
    temp4$Values[if(length(m)>0){-m}else{1:5}]<-rvs$bar[na.omit(match(col_anno[,1],names(rvs$bar)))]
    #####barplot
    rvs$counter->l
    ggplot(temp4, aes(y=Values, x=State)) + 
      geom_bar(stat="identity",fill=temp4$Col) +
      scale_y_continuous(trans = 'log10') +
      ggtitle(paste("Cell counts after",l,"hour(s)"))+
      labs(x='State',y='number of cells')
  })
  observe({
    mt_l<-input$mt_l
    mt_c<-input$mt_c
    #####
    ret3_f<-function(mt_l,mt_c){
      mat.or.vec(mt_l,mt_c)->y
      y[1:length(y)]<-c(1:length(y))
      m2<-cbind(NA,rbind(NA,y,NA),NA)
      addresses <- expand.grid(x = 1:nrow(y), y = 1:ncol(y))
      ret<-c()
      for(i in 1:-1)
        for(j in 1:-1)
          if(i!=0 || j !=0)
            ret<-rbind(ret,m2[addresses$x+i+1+nrow(m2)*(addresses$y+j)])
      y2<-cbind(y[,ncol(y)],y)
      y2<-cbind(y2,y[,1])
      m3<-cbind(NA,rbind(NA,y2,NA),NA)
      addresses2 <- expand.grid(x = 1:nrow(y2), y = 1:ncol(y2))
      ret2<-c()
      for(i in 1:-1)
        for(j in 1:-1)
          if(i!=0 || j !=0)
            ret2<-rbind(ret2,m3[addresses2$x+i+1+nrow(m3)*(addresses2$y+j)])
      ret2[,-c(1:nrow(y2),(1+(ncol(y2)-1)*nrow(y2)):(ncol(y2)*nrow(y2)))]->ret3
      return(ret3)
    }
    rvs$ret3<-ret3_f(mt_l,mt_c)
    rvs$iter=array("S",dim=c(mt_l,mt_c))
    rvs$bar=setNames(c(mt_l*mt_c,0,0,0,0),c("S","E","I","R","D"))
  })
  #########define the update observations
  observe({
    ############asign variables based on user input
    vD<-input$vD
    d0<-input$d0
    vT<-input$vT
    TD<-input$TD
    death_rate<-input$death_rate
    diffusion<-input$diffusion
    n_iter<-input$n_iter
    mt_l<-input$mt_l
    mt_c<-input$mt_c
    ############define CA update rules based on the user input and SEIR model
    state_evolve<-function(state){
      if(state=="S"){
        if(rexp(1,vD)<1){state<-"E"}
        else{state="S"}
      }
      else
        if(state=="E"){
          if(rexp(1,TD*vT+d0)<1){
            if(rexp(1,TD*vT)<rexp(1,d0)){state="I"}else{state="S"}
          }
          else{state="E"}
        }
      else
        if(state=="I"){
          if(rexp(1,death_rate+d0)<1){
            if(rexp(1,death_rate)<rexp(1,d0)){state="D"}else{state="R"}
          }
          else{state="I"}
        }
      else
        if(state=="R"){
          if(rexp(1,death_rate+vD)<1){
            if(rexp(1,death_rate)<rexp(1,vD)){state="D"}else{state="E"}
          }
          else{state="R"}
        }
      return(state)
    }
    #############Diffusion interaction model on neighbours
    diffusion_function<-function(state){
      if(state=="S"){
        if(rexp(1,diffusion)<1){state<-"R"}
        else{state="S"}
      }
      if(state=="E"){
        if(rexp(1,diffusion)<1){state<-"I"}
        else{state="E"}
      }
      return(state)
    }
    ###################
    #################if statement to hold the CA implimentation until we click Go, gives user time to select params
    if(rvs$reset==0){
      return(NULL)
    }
    #################initiate model
    isolate({
      ######
      ######
      if(rvs$counter==0){
        rvs$iter<-array("S",dim=c(mt_l,mt_c))
        rvs$bar<-setNames(c(mt_l*mt_c,0,0,0,0),c("S","E","I","R","D"))
      }
      ########define update function, counter is dynamic variable which counts the number of hours we implement
      if(rvs$counter>0){
        rvs$iter->temp2
        ##run intracellular dynamic defined by state update function
        mapply(function(j){return(mapply(function(i){return(state_evolve(temp2[i,j]))},1:nrow(temp2)))},1:mt_c)->temp2
        ###run neighbour interaction update based on diffusion function
        c(which(temp2=="I"),which(temp2=="R"))->temp
        if(length(temp)>0){
          for(k in 1:length(temp)){
            temp2[as.vector(na.omit(rvs$ret3[,temp[k]]))]->tt
            mapply(function(i)return(diffusion_function(tt[i])),1:length(tt))->temp2[as.vector(na.omit(rvs$ret3[,temp[k]]))]
          }
        }
        ####asign updated matrix to dynamic
        temp2->rvs$iter
        table(temp2)->rvs$bar
      }
      ####update number of hours run
      rvs$counter<-rvs$counter+1
    })
    ####stop when we reach n_iter hours - user input for lenght of simulation
    if (((isolate(rvs$counter<n_iter))))
      invalidateLater(200,session)
  })
  #################allow action button to start automaton
  observe({
    if(input$go>0){
      rvs$reset<<-1
    }
    
  })
}