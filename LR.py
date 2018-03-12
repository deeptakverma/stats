import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

feature_data = pd.read_csv("data.csv")
print(feature_data)

X = feature_data.drop('Activity', axis=1)
print(X)

lm = LinearRegression()
lm.fit(X,feature_data.Activity)

print("intercept", lm.intercept_)
print("# of coeff", len(lm.coef_))
#print(X.columns[1])
#print(lm.coef_[1])
for i in range(len(X.columns)): print(X.columns[i], lm.coef_[i])

#pd.DataFrame(zip(X.columns, lm.coef_), columns = ['features','estCoeffs'])
