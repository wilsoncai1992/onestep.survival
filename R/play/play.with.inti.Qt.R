# examine initial estimate for conditional survival curve
play <- colMeans(Qn.A1.t)
# play <- play *10
lines(play~ T.uniq, col = 'green')
