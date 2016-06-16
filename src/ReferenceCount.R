public = function(x) {
  returnVal = data.frame("ReferenceCount" = integer(0));
  tmp = data.frame("Run" = character(0), stringsAsFactors=FALSE);
  count = 1;
  index = 1;
  for (i in x[[1]]) {
    print(count)
    if (inherits(try(getReference(i), silent = TRUE),
                 "try-error") || dim(getReference(i))[1] < 1) {
    } else {
      tmp[index, ] = c(i)
      returnVal[index, ] = dim(getReference(i))[1]
      index = index + 1
    }
    count = count + 1
  }
  returnVal = cbind(tmp, returnVal)
  returnVal
}