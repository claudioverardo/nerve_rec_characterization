
def get_fiber_model_id(fiber_model):

    if isinstance(fiber_model, str):

        if fiber_model == "ascent_sundt15":
            fiber_model_id = 1
        elif fiber_model == "ascent_mrg":
            fiber_model_id = 6
        else:
            raise Exception("fiber_model not valid")

    else:

        if fiber_model == 1:
            fiber_model_id = "ascent_sundt15"
        elif fiber_model == 6:
            fiber_model_id = "ascent_mrg"
        else:
            raise Exception("fiber_model_id not valid")

    return fiber_model_id


def check_fiber_myelinated(fiber_model):
    
    return fiber_model == "ascent_mrg"


def check_fiber_unmyelinated(fiber_model):
    
    return fiber_model == "ascent_sundt15"
