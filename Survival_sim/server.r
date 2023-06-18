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

load("files_for_KM_plot_shiny_not_refined.rd")
library(shiny)
library("tidyverse")
library('deSolve')
library('survival')
library("survminer")
library('ggplot2')
server<-function(input,output,session){
  #########define reative values which update as the automaton evolves
  rvs<-reactiveValues(alive=array(c(0,0,0),dim=c(100*24,3)),D4n=array(c(0,0,0),dim=c(100*24,3)),
                      D4Tn=array(c(0,0,0),dim=c(100*24,3)),events_a=array(c(0,0,0),dim=c(20000,3)),
                      events_D4n=array(c(0,0,0),dim=c(20000,3)),events_D4Tn=array(c(0,0,0),dim=c(20000,3)),
                      reset=0)
  #########output dynamically updating plot of automaton at the top
  #########define the update observations
  observe({
    ############asign variables based on user input
    vD<-input$vD
    d0<-input$d0
    vT<-input$vT
    TD<-input$TD
    death_rate<-input$death_rate
    #define sc_model
    sc_model <- function (t, x, params) {
      ## state variables
      S <- x[1]
      E <- x[2]
      I <- x[3]
      R <- x[4]
      D <- x[5]
      ## parameters
      vD <- params["vD"]
      d0 <- params["d0"]
      vT <- params["vT"]
      TD <- params["TD"]
      death_rate <- params["death_rate"]
      N <- S+E+I+R
      ## now code the model equations
      dSdt <- d0*E - vD*S
      dEdt <- vD*S - (d0+TD*vT)*E
      dIdt <- TD*vT*E + vD*R - d0*I - death_rate*I
      dRdt <- d0*I - vD*R - death_rate*R
      dDdt <- death_rate*(I+R)
      ## combine results into a single vector
      dxdt <- c(dSdt,dEdt,dIdt,dRdt,dDdt)
      ## return result as a list
      list(dxdt)
    }
    ###############execute for true params
    parms <- c(d0=0.2464819,vD=0.002113848,vT=6.41248,death_rate=1/20.2,TD=1/13)
    times <- seq(from=0,to=24*100,by=1)
    xstart <- c(S=5133*(1+per_D_at_D3_sc_adj_trans),E=0,I=0,R=0,D=0)
    ######
    ode(
      func=sc_model,
      y=xstart,
      times=times,
      parms=parms
    ) %>%
      as.data.frame() -> out
    ######convert this to %alive at each time point
    5133*(1+per_D_at_D3_sc_adj_trans)->n_sc
    1-out$D/n_sc->per_alive_true
    ######now to %DUX4 naieve at each time point
    out$S/n_sc->per_D4n_true
    ######now to %target naieve at each time point
    (out$S+out$E)/n_sc->per_D4Tn_true
    ##################now run for a perturbed parameter
    if(d0==0.2464819 & vD==0.002113848 & vT==6.41248 & death_rate==1/20.2 & TD==1/13){out_pert<-out}else{parms <- c(d0=d0,vD=vD,vT=vT,death_rate=death_rate,TD=TD)
    ###############
    ode(
      func=sc_model,
      y=xstart,
      times=times,
      parms=parms
    ) %>%
      as.data.frame() -> out_pert}
    ############
    1-out_pert$D/n_sc->per_alive_p
    out_pert$S/n_sc->per_D4n_p
    (out_pert$S+out_pert$E)/n_sc->per_D4Tn_p
    data.frame(out$time,per_alive_true,per_alive_p)->rvs$alive
    data.frame(out$time,per_D4n_true,per_D4n_p)->rvs$D4n
    data.frame(out$time,per_D4Tn_true,per_D4Tn_p)->rvs$D4Tn
  })
  output$plot1<-renderPlot({
    rvs$alive->alive_comp
    names(alive_comp)<-c("time","True Parameters","Test Parameters")
    gather(alive_comp,variable,value,-time)->alive
    #######
    ggplot(alive,aes(x=time,y=value,color=variable))+
      geom_line(size=2)+
      theme_classic()+
      labs(x='time (hrs)',y='% cells alive')+
      ggtitle("Proportion of Cells Alive Single cell model")
  })
  #####automaton
  observe({
    if(rvs$reset==0){
      return(NULL)
    }
    ####start with automaton here
    ###define inputs for true model first
    d0=0.2464819
    vD=0.002113848
    vT=6.41248
    death_rate=1/20.2
    TD=1/13
    diffusion=0.0402
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
    ############
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
    ret3<-ret3_f(mt_l,mt_c)
    ############
    ####initialise
    array("S",dim=c(mt_l,mt_c))->iter
    CA_D4Tn_true<-CA_D4n_true<-CA_survival_true<-array(c(0,0),dim=c(length(iter),2))
    ###########################loop iterate the remaining proportions
    for(i in 1:n_iter){
      ##run intracellular dynamic
      mapply(function(j){return(mapply(function(i){return(state_evolve(iter[i,j]))},1:mt_l))},1:mt_c)->iter
      ###run interaction update
      c(which(iter=="I"),which(iter=="R"))->temp
      if(length(temp)>0){
        for(k in 1:length(temp)){
          iter[as.vector(na.omit(ret3[,temp[k]]))]->tt
          mapply(function(i)return(diffusion_function(tt[i])),1:length(tt))->iter[as.vector(na.omit(ret3[,temp[k]]))]
        }
      }
      ######output survival data
      which(CA_survival_true[,1]==1)->u
      which(iter=="D")->v
      if(length(u>0)){match(u,v)->w;v[-w]->v}
      CA_survival_true[v,2]<-i
      CA_survival_true[v,1]<-1
      #######
      which(CA_D4n_true[,1]==1)->u
      which(iter!="S")->v
      if(length(u>0)){match(u,v)->w;v[-na.omit(w)]->v}
      CA_D4n_true[v,2]<-i
      CA_D4n_true[v,1]<-1
      # #######
      which(CA_D4Tn_true[,1]==1)->u
      which(iter!="S" & iter!="E")->v
      if(length(u>0)){match(u,v)->w;v[-na.omit(w)]->v}
      CA_D4Tn_true[v,2]<-i
      CA_D4Tn_true[v,1]<-1
    }
    ####now consider the perturbed params
    ############
    ############
    vD<-input$vD
    d0<-input$d0
    vT<-input$vT
    TD<-input$TD
    death_rate<-input$death_rate
    diffusion<-input$diffusion
    mt_l<-input$mt_l
    mt_c<-input$mt_c
    ############
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
    ret3<-ret3_f(mt_l,mt_c)
    ############
    ####initialise
    array("S",dim=c(mt_l,mt_c))->iter
    ############
    ######define matrices to keep track of the fate of cells
    CA_D4Tn_p<-CA_D4n_p<-CA_survival_p<-array(c(0,0),dim=c(length(iter),2))
    ###########################loop iterate the remaining proportions
    for(i in 1:n_iter){
      ##run intracellular dynamic
      mapply(function(j){return(mapply(function(i){return(state_evolve(iter[i,j]))},1:mt_l))},1:mt_c)->iter
      ###run interaction update
      c(which(iter=="I"),which(iter=="R"))->temp
      if(length(temp)>0){
        for(k in 1:length(temp)){
          iter[as.vector(na.omit(ret3[,temp[k]]))]->tt
          mapply(function(i)return(diffusion_function(tt[i])),1:length(tt))->iter[as.vector(na.omit(ret3[,temp[k]]))]
        }
      }
      ######plot automata
      which(CA_survival_p[,1]==1)->u
      which(iter=="D")->v
      if(length(u>0)){match(u,v)->w;v[-w]->v}
      CA_survival_p[v,2]<-i
      CA_survival_p[v,1]<-1
      #################
      which(CA_D4n_p[,1]==1)->u
      which(iter!="S")->v
      if(length(u>0)){match(u,v)->w;v[-na.omit(w)]->v}
      CA_D4n_p[v,2]<-i
      CA_D4n_p[v,1]<-1
      # #######
      which(CA_D4Tn_p[,1]==1)->u
      which(iter!="S" & iter!="E")->v
      if(length(u>0)){match(u,v)->w;v[-na.omit(w)]->v}
      CA_D4Tn_p[v,2]<-i
      CA_D4Tn_p[v,1]<-1
    }
    library('survival')
    library("survminer")
    ####################
    ####################
    #######
    c(rep("1",nrow(CA_survival_true)),rep("0",nrow(CA_survival_p)))->l
    rbind(CA_survival_true,CA_survival_p)->events
    data.frame(events)->events
    names(events)<-c("event","time")
    events$model<-l
    events->rvs$events_a
    ###############
    c(rep("1",nrow(CA_survival_true)),rep("0",nrow(CA_survival_p)))->l
    rbind(CA_D4n_true,CA_D4n_p)->events
    data.frame(events)->events
    names(events)<-c("event","time")
    events$model<-l
    events->rvs$events_D4n
    # ###################
    c(rep("1",nrow(CA_survival_true)),rep("0",nrow(CA_survival_p)))->l
    rbind(CA_D4Tn_true,CA_D4Tn_p)->events
    data.frame(events)->events
    names(events)<-c("event","time")
    events$model<-l
    events->rvs$events_D4Tn
  })
  
  observe({
    if(input$go>0){
      rvs$reset<<-1
    }
  })
  ###############
  output$plot4<-renderPlot({
    if(rvs$reset==0){
      return(NULL)
    }
    rvs$events_a->events
    survfit(Surv(time,event)~model,data=events)->fit
    #hist(events$time)
    ggsurvplot(fit,data=events,
               conf.int = TRUE,          # Add confidence interval
               pval = TRUE,              # Add p-value
               legend.labs =
                 c("Test Parameters", "True Parameters"),
               title="Proportion of Cells Alive Cell Automaton model")
  })
}
