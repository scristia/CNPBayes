test_example <- function(){
  myfunction <- function(x, y){
    x+y
  }
  result <- myfunction(2, 2)
  checkIdentical(result, 4)
  result <- myfunction(2.1, pi)
  checkEquals(result, 2.1+3.142, tolerance=0.005)
  checkException(myfunction("a", "b"))
}
