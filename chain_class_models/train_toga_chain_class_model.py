#!/usr/bin/env python3
"""Script to train XGBoost models."""
import os
from datetime import datetime as dt
import xgboost as xgb
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
import joblib

# lists of features required by single and multi exon models
SE_MODEL_FEATURES = ["gl_exo", "flank_cov", "exon_perc", "synt_log"]
ME_MODEL_FEATURES = ["gl_exo", "loc_exo", "flank_cov", "synt_log", "intr_perc"]


def train_on(x, y, save_to, name=None):
    """Train model on the X and y given."""
    # models parameters, work fine for both multi and single exon models
    n_trees = 50
    max_depth = 3
    learning_rate = 0.1
    # create and fit the model, also add cross-validation
    model = xgb.XGBClassifier(
        n_estimators=n_trees, max_depth=max_depth, learning_rate=learning_rate
    )
    kfold = StratifiedKFold(n_splits=5, random_state=777, shuffle=True)
    results = cross_val_score(model, x, y, cv=kfold)
    model.fit(x, y)
    if name:  # some verbosity
        y_lst = list(y)
        print(f"{name} model: ")
        print(f"Training on {len(x)} samples")
        print(f"Positives: {y_lst.count(1)}; Negatives: {y_lst.count(0)}")
        print(f"Using features: {x.columns}")
    print("Accuracy: {0:.3f} {1:.3f}".format(results.mean() * 100, results.std() * 100))
    joblib.dump(model, save_to)  # save the model
    print(f"Model saved to: {save_to}")


# load dataset, defile where is what: input and output files
# training data is in the repository
t0 = dt.now()
file_location = os.path.dirname(__file__)
models_dir = "."
train_tsv = os.path.join(file_location, "train.tsv")
se_model_dat = os.path.join(file_location, "se_model.dat")
me_model_dat = os.path.join(file_location, "me_model.dat")

# load dataframe
df = pd.read_csv(train_tsv, header=0, sep="\t")
print(f"Training dataset size: {len(df)}")

# split overall dataframe into two frames:
# one for 1-exon genes, another for multi-exon
df_se = df[df["single_exon"] == 1]
df_me = df[df["single_exon"] == 0]

# create X and y for both models
# Single and multi-exon models require different sets of features
X_se = df_se.copy()
X_se = X_se[SE_MODEL_FEATURES]
y_se = df_se["y"]

X_me = df_me.copy()
X_me = X_me[ME_MODEL_FEATURES]
y_me = df_me["y"]

print(f"Single exon train length: {len(X_se)}; multi exon: {len(X_me)}")

# train and save models
train_on(X_se, y_se, se_model_dat, name="Single exon")
train_on(X_me, y_me, me_model_dat, name="Multi exon")

print(f"Done in {dt.now() - t0}")
