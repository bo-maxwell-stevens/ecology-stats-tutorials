# Random Forest

`random_forest.py` uses a group-aware 80/20 split and fits imputation only on training data. R-squared compares predictions with a mean baseline, so it can be negative; RMSE and MAE are errors in response units. Importance is descriptive, not causal.

`leave_one_site_out.py` uses `LeaveOneGroupOut` and reports **pooled** metrics, never invalid one-observation site R-squared values. In the supplied metadata every site has one record, so LOSO is effectively leave-one-record-out; use field/site groups with multiple observations for a genuine unobserved-site test. Run both scripts with `python 07_Random_Forest/<script>.py`. Documentation: [scikit-learn](https://scikit-learn.org/stable/).
