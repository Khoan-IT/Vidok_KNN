from sklearn.neighbors import NearestNeighbors

def create_model(model_name, *args, **kwargs):
    model = None

    if model_name == "nearest_neighbors":
        model = NearestNeighbors(*args, **kwargs)

    return model
    