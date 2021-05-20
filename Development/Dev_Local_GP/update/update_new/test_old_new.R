load(file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/joint_model_fit_1_old.RData'))
load(file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/joint_model_fit_1_new.RData'))

all.equal(joint_model_fit_1_new, joint_model_fit_1_old)

joint_model_fit_1_new$running_time
joint_model_fit_1_old$running_time
