
from collagen_amp_ml.train import train_amp_classifier, fit_jnp_domain_adapter
print(train_amp_classifier("data/public/veltri_positive.csv","data/public/veltri_negative.csv","models"))
print(fit_jnp_domain_adapter("data/jnp_collagen_amp_ml_anchors.csv","models"))
