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
library('ggplot2')
load("files_for_SEIR_shiny_not_refined.rd")
server<-function(input,output,session){
  rvs<-reactiveValues(out_sc=NULL,out_sn=NULL)
  observe({
    vD<-input$vD
    d0<-input$d0
    vT<-input$vT
    TD<-input$TD
    death_rate<-input$death_rate
    diffusion<-input$diffusion
    n_iter<-input$n_iter
    n_sc<-input$n_sc
    n_sn<-input$n_sn
    ###############
    ###############
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
    parms <- c(d0=d0,vD=vD,vT=vT,death_rate=death_rate,TD=TD)
    times <- seq(from=0,to=n_iter,by=1)
    xstart <- c(S=n_sc,E=0,I=0,R=0,D=0)
    ###
    library(tidyverse)
    ode(
      func=sc_model,
      y=xstart,
      times=times,
      parms=parms
    ) %>%
      as.data.frame() -> rvs$out_sc
    ###########SN model
    sn_model <- function (t, x, params) {
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
      diffusion <- params["diffusion"]
      N <- S+E+I+R
      ## now code the model equations
      dSdt <- d0*E - vD*S - diffusion*(I+R)*S
      dEdt <- vD*S - (d0+TD*vT)*E - diffusion*(I+R)*E
      dIdt <- TD*vT*E + vD*R - d0*I - death_rate*I + diffusion*(I+R)*E
      dRdt <- d0*I - vD*R - death_rate*R + diffusion*(I+R)*S
      dDdt <- death_rate*(I+R)
      ## combine results into a single vector
      dxdt <- c(dSdt,dEdt,dIdt,dRdt,dDdt)
      ## return result as a list!
      list(dxdt)
    }
    parms2 <- c(d0=d0,vD=vD,vT=vT,death_rate=death_rate,diffusion=diffusion,TD=TD)
    times2 <- seq(from=0,to=n_iter,by=1)
    xstart2 <- c(S=n_sn,E=0,I=0,R=0,D=0)
    ###
    ode(
      func=sn_model,
      y=xstart2,
      times=times2,
      parms=parms2
    ) %>%
      as.data.frame() -> rvs$out_sn
  })
  output$plot1<-renderPlot({
    rvs$out_sc -> out
    gather(out,variable,value,-time)->temp5
    temp5$variable<-factor(temp5$variable,levels=col_anno[,1])
    ggplot(temp5,aes(x=time,y=value,color=variable))+
      geom_line(size=2)+
      scale_color_manual(values=col_anno[,2])+
      labs(x='time (hrs)',y='number of cells')+
      ggtitle("Temporal Evolution: FSHD myocytes")
  })
  output$plot2<-renderPlot({
    rvs$out_sc->out
    data.frame(rep(c("S","E","I","R"),2),c(rep("Model",4),rep("Data",4)),
               c(as.numeric(out[which(out$time==24*3),2:5]),4956,14,13,150))->df2
    names(df2)<-c("State","Method","Values")
    ###
    ggplot(df2, aes(fill=Method, y=Values, x=State)) + 
      geom_bar(position="dodge", stat="identity") +
      scale_y_continuous(trans = 'log10')+
      labs(x='State',y='number of cells')+
      ggtitle("Comparison of Day 3 of model with scRNAseq of FSHD myocytes")
  })
  ############
  ############
  output$plot3<-renderPlot({
    rvs$out_sn -> out2
    gather(out2,variable,value,-time)->temp6
    temp6$variable<-factor(temp6$variable,levels=col_anno[,1])
    ggplot(temp6,aes(x=time,y=value,color=variable))+
      geom_line(size=2)+
      scale_color_manual(values=col_anno[,2])+
      labs(x='time (hrs)',y='number of nuclei')+
      ggtitle("Temporal Evolution: FSHD myotubes")
  })
  output$plot4<-renderPlot({
    rvs$out_sn->out2
    data.frame(rep(c("S","E","I","R"),2),c(rep("Model",4),rep("Data",4)),
               c(as.numeric(out2[which(out2$time==24*2),2:5]),58,0,3,78))->df2
    names(df2)<-c("State","Method","Values")
    ###
    ggplot(df2, aes(fill=Method, y=Values, x=State)) + 
      geom_bar(position="dodge", stat="identity") +
      #scale_y_continuous(trans = 'log10')+
      labs(x='State',y='number of cells')+
      ggtitle("Comparison of Day 2 of model with snRNAseq of FSHD myocytes")
  })
}
