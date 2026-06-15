#!/bin/bash
# Simple AIC and DeltaL calculator
# Usage: 
#   aic.sh k lhood              (AIC only)
#   aic.sh k est obs            (AIC + DeltaL)

k=$1
est=$2
obs=$3

if [ $# -eq 2 ]; then
    # 只计算 AIC
    awk -v k="$k" -v est="$est" 'BEGIN {
        ln10 = log(10);
        aic = 2*k - 2*(est * ln10);
        printf "AIC = %.4f\n", aic;
        printf "ln(L) = %.6f\n", est * ln10;
    }'
elif [ $# -eq 3 ]; then
    # 计算 AIC 和 DeltaL
    awk -v k="$k" -v est="$est" -v obs="$obs" 'BEGIN {
        ln10 = log(10);
        aic = 2*k - 2*(est * ln10);
        deltaL = obs - est;
        printf "AIC = %.4f\n", aic;
        printf "ln(L_est) = %.6f\n", est * ln10;
        printf "DeltaL = %.4f\n", deltaL;
        if (deltaL > 50) printf "Warning: Large DeltaL!\n";
    }'
else
    echo "Usage: aic.sh k est [obs]"
    exit 1
fi
