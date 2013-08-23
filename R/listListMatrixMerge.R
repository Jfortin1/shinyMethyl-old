listListMatrixMerge <-
function(listList1,listList2){
		newList <- vector("list",length(listList2))
		names(newList) <- names(listList2)
		for (k in 1:length(listList2)){
			newList[[k]] <- listMatrixMerge(listList1[[k]],listList2[[k]])
		}
		return(newList)
	}
