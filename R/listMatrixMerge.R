listMatrixMerge <-
function(list1,list2){
		newList <- vector("list",length(list2))
		names(newList) <- names(list2)
		for (k in 1:length(list2)){
			newList[[k]] <- cbind(list1[[k]],list2[[k]])
		}
		return(newList)
	}
