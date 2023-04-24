from django.shortcuts import render
from django.conf import settings
from .forms import ParameterForm
# Create your views here.
def index(request):

    context ={}
    context['form']= ParameterForm()

    return render(request, str(settings.BASE_DIR) + "/static/templates/index.html", context)