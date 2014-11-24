#!/usr/bin/env Rscript

n=commandArgs(trailing=TRUE)
##Yes, this should all be redone with tidy/dplyr, etc.,
##But I'm trying not to introduce any new dependencies.
#output=read.table(n[1],header=T,colClasses=c("character",rep("integer",10)))
output = read.table(n[1],header=T,stringsAsFactors=FALSE)
native=read.table(n[2],stringsAsFactors=FALSE)#,colClasses=c("character",rep("integer",2)))
truth=read.table(n[3],colClasses=c("character",rep("integer",2)),stringsAsFactors=FALSE)
outfilename=n[4]

renameOutput = function(x)
    {
        names(x) <- c("chromo","pfirst","plast","chrom2","mfirst","mlast","name","score","strand1","strand2","nplus","nminus","pdist","pin","mdist","min")
        x
    }

output = renameOutput(output)
FOUNDNATIVE=array(NA,dim=nrow(output))
FOUNDSIM=array(NA,dim=nrow(output))
for( O in 1:nrow(output) )
    {
        for(n in 1:nrow(native))
            {
                if ( output$chromo[O] == native$V1[n] ) ##Same chromo?
                    {
                        if( output$pfirst[O] != -1 )
                            {
                                if( (output$pfirst[O] >= native$V2[n] &
                                     output$pfirst[O] <= native$V3[n] ) |
                                   (output$plast[O] >= native$V2[n] &
                                    output$plast[O] <= native$V3[n] ) )
                                    {
                                        FOUNDNATIVE[O]=n
                                    }
                            } else if ( output$mfirst[O] != -1 )
                                {
                                    if( (output$mfirst[O] >= native$V2[n] &
                                         output$mfirst[O] <= native$V3[n] ) |
                                       (output$mlast[O] >= native$V2[n] &
                                        output$mlast[O] <= native$V3[n] ) )
                                        {
                                            FOUNDNATIVE[O]=n
                                        }
                                }
                    }
            }
        for(n in 1:nrow(truth))
            {
                if ( output$chromo[O] == truth$V1[n] ) ##Same chromo?
                    {
                        if( output$pfirst[O] != -1 )
                            {
                                A = truth$V2[n] - output$pfirst[O]
                                B = truth$V2[n] - output$plast[O] 
                                if( (A>0 & A<= 250) |
                                   (B>0 & B<= 250 ) )
                                    {
                                        FOUNDSIM[O]=n
                                    }
                            }
                        if ( output$mfirst[O] != -1 )
                            {
                                A = output$mfirst[O] - truth$V2[n]
                                B = output$mlast[O] - truth$V2[n]
                                if( (A>0 & A<= 250) |
                                   (B>0 & B<= 250 ) )
                                    {
                                        FOUNDSIM[O]=n
                                    }
                            }
                    }
            }
    }

OO=as.data.frame(cbind(FOUNDNATIVE,FOUNDSIM,output))
write.table(OO,file=outfilename,row.names=F,quote=F)
