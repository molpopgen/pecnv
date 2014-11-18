##Yes, this should all be redone with tidy/dplyr, etc.,
##But I'm trying not to introduce any new dependencies.
output=read.table("teclust_output.gz",header=T,colClasses=c("character",rep("integer",10)))
native=read.table("TE_position_r5.1",colClasses=c("character",rep("integer",2)))
truth=read.table("line99_truth",colClasses=c("character",rep("integer",2)))

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
                                        print(paste("plus",truth$V2[n],output$pfirst[O],output$plast[O]))
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
                                        print(paste("minus",truth$V2[n],output$mfirst[O],output$mlast[O]))
                                        FOUNDSIM[O]=n
                                    }
                            }
                    }
            }
    }

